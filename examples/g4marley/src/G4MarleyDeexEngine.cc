/// @file
/// @brief Thin wrapper around marley::NucleusDecayer for Geant4 fragments
//
// This file is part of the MARLEY Geant4 de-excitation interface (g4marley).
// It is distributed under the same terms as MARLEY (GPL v3).

// Standard library includes
#include <chrono>
#include <map>
#include <cmath>
#include <random>
#include <stdexcept>

// Geant4 includes
#include "G4Fragment.hh"
#include "G4Neutron.hh"
#include "G4NucleiProperties.hh"
#include "G4Proton.hh"
#include "G4Electron.hh"
#include "G4EventManager.hh"
#include "G4HadronicProcess.hh"
#include "G4Nucleus.hh"
#include "G4SteppingManager.hh"
#include "G4Track.hh"
#include "G4TrackingManager.hh"
#include "G4ios.hh"
#include "Randomize.hh"

// HepMC3 includes
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

// MARLEY includes
#include "marley/Error.hh"
#include "marley/Generator.hh"
#include "marley/JSON.hh"
#include "marley/JSONConfig.hh"
#include "marley/LevelDensityModel.hh"
#include "marley/OpticalModel.hh"
#include "marley/MassTable.hh"
#include "marley/NucleusDecayer.hh"
#include "marley/Parity.hh"
#include "marley/StructureDatabase.hh"
#include "marley/hepmc3_utils.hh"
#include "marley/marley_utils.hh"

#include "G4MarleyDeexEngine.hh"

namespace {

  // Speed of light in cm/ns. MARLEY stores decay times as c*t in cm.
  constexpr double C_CM_PER_NS = 29.9792458;

  // Tolerance (MeV) for the energy-conserving momentum rescaling
  constexpr double E_TOLERANCE = 1e-9;

  // Reads an integer attribute, returning a default if it is absent
  int int_attr( const HepMC3::ConstGenParticlePtr& p, const std::string& name,
    int def )
  {
    auto a = p->attribute< HepMC3::IntAttribute >( name );
    return a ? a->value() : def;
  }

  double double_attr( const HepMC3::ConstGenParticlePtr& p,
    const std::string& name, double def )
  {
    auto a = p->attribute< HepMC3::DoubleAttribute >( name );
    return a ? a->value() : def;
  }

  // Geant4 mass (MeV) for a particle MARLEY may emit. Returns a negative
  // value for anything unexpected.
  double g4_mass( int pdg ) {
    if ( pdg == marley_utils::PHOTON ) return 0.;
    if ( pdg == marley_utils::NEUTRON ) return G4Neutron::Neutron()->GetPDGMass();
    if ( pdg == marley_utils::PROTON ) return G4Proton::Proton()->GetPDGMass();
    if ( pdg == marley_utils::ELECTRON ) {
      return G4Electron::Electron()->GetPDGMass();
    }
    if ( marley_utils::is_ion(pdg) ) {
      int Z = marley_utils::get_particle_Z( pdg );
      int A = marley_utils::get_particle_A( pdg );
      return G4NucleiProperties::GetNuclearMass( A, Z );
    }
    return -1.;
  }

}

G4MarleyDeexEngine::G4MarleyDeexEngine( const G4MarleyDeexConfig& config )
  : config_( config )
{
  // A JSON config with null "reactions" and "source" makes MARLEY build only
  // what de-excitation needs (mass table, structure database, RNG). This is
  // the same pattern used by "marley decay" and marley_deexcite_hepmc.
  marley::JSON json = marley::JSON::object();
  json["reactions"] = marley::JSON( nullptr );
  json["source"] = marley::JSON( nullptr );

  std::uint64_t seed = config_.seed;
  if ( seed == 0 ) {
    // Derive MARLEY's seed from the Geant4 engine so that a whole run is
    // reproducible from the Geant4 seed. Keep it positive as a signed long.
    std::uint64_t hi = static_cast< std::uint64_t >( G4UniformRand() * 2147483647. );
    std::uint64_t lo = static_cast< std::uint64_t >( G4UniformRand() * 4294967295. );
    seed = ( hi << 32 ) | lo;
    if ( seed == 0 ) seed = 1;
  }
  config_.seed = seed;
  json["seed"] = static_cast< long long >( seed );

  marley::JSONConfig jc( json );
  gen_ = std::make_unique< marley::Generator >( jc.create_generator() );

  G4cout << "[g4marley] MARLEY de-excitation engine ready (seed " << seed
    << ", spin model " << static_cast<int>( config_.spin_model )
    << ", min Ex " << config_.min_Ex << " MeV)" << G4endl;
}

G4MarleyDeexEngine::~G4MarleyDeexEngine() = default;

G4MarleyDeexEngine::Skip G4MarleyDeexEngine::check( const G4Fragment& frag )
  const
{
  if ( !config_.use_marley ) return Skip::Disabled;
  const int A = frag.GetA_asInt();
  const double Ex = frag.GetExcitationEnergy();
  if ( A < config_.min_A ) return Skip::Light;
  if ( Ex < config_.min_Ex ) return Skip::LowEx;
  if ( Ex > config_.max_Ex_per_A * A ) return Skip::HighEx;
  return Skip::None;
}

bool G4MarleyDeexEngine::decay_with_attempts( const G4Fragment& frag,
  G4MarleyDeexRecord& rec, std::string& last_error )
{
  auto& recorder = G4MarleyDeexRecorder::Instance();
  for ( int attempt = 1; attempt <= config_.max_attempts; ++attempt ) {
    rec = G4MarleyDeexRecord();
    rec.Z = frag.GetZ_asInt();
    rec.A = frag.GetA_asInt();
    rec.Ex = frag.GetExcitationEnergy();
    rec.p4 = frag.GetMomentum();
    rec.attempts = attempt;
    bool ok = false;
    try {
      ok = this->try_decay_once( frag, rec );
      if ( !ok ) last_error = "energy/particle bookkeeping failed";
    }
    catch ( const std::exception& e ) {
      // marley::Error derives from std::runtime_error
      last_error = e.what();
    }
    if ( ok ) return true;
    ++recorder.counters().failed_attempts;
  }
  return false;
}

bool G4MarleyDeexEngine::unwanted_cn_outcome( const G4MarleyDeexRecord& rec )
  const
{
  if ( rec.spin_source != 2 || rec.products.empty() ) return false;
  const auto& residual = rec.products.back();
  const size_t n_emitted = rec.products.size() - 1;

  // Compound elastic: only the projectile comes out, target left in its
  // ground state (MARLEY always ends the cascade in the ground state)
  if ( config_.cn_reject_elastic && n_emitted == 1
    && rec.products.front().pdg == rec.proj_pdg
    && residual.pdg == rec.target_pdg ) return true;

  // Radiative capture of a neutron: compound nucleus survives, only gammas
  if ( config_.cn_reject_neutron_capture
    && rec.proj_pdg == marley_utils::NEUTRON
    && residual.pdg == marley_utils::get_nucleus_pid(rec.Z, rec.A) )
  {
    bool only_gammas = true;
    for ( size_t i = 0; i < n_emitted; ++i ) {
      int pdg = rec.products[i].pdg;
      if ( pdg != marley_utils::PHOTON && pdg != marley_utils::ELECTRON ) {
        only_gammas = false;
      }
    }
    if ( only_gammas ) return true;
  }
  return false;
}

bool G4MarleyDeexEngine::decay( const G4Fragment& frag,
  G4MarleyDeexRecord& rec )
{
  auto& counters = G4MarleyDeexRecorder::Instance().counters();
  auto start = std::chrono::steady_clock::now();

  std::string last_error;
  bool ok = this->decay_with_attempts( frag, rec, last_error );

  // Condition compound-nucleus decays on the channels the triggering
  // cross section actually contains. Rejection sampling over (J, decay)
  // jointly gives P(J, decay | accepted), as it should.
  int n_resample = 0;
  if ( ok && this->unwanted_cn_outcome(rec) ) {
    const auto key = std::make_tuple( rec.target_pdg, rec.proj_pdg,
      std::lround(rec.proj_KE * 1000.) );
    const bool known_closed = ( closed_channels_.count(key) > 0 );
    const int max_resample = known_closed ? 10 : config_.cn_max_resample;
    while ( this->unwanted_cn_outcome(rec) && n_resample < max_resample ) {
      G4MarleyDeexRecord trial;
      ++n_resample;
      if ( this->decay_with_attempts(frag, trial, last_error) ) rec = trial;
    }
    if ( this->unwanted_cn_outcome(rec) ) {
      ++counters.cn_resample_failed;
      closed_channels_.insert( key );
    }
  }
  rec.cn_resamples = n_resample;
  counters.cn_resamples += n_resample;

  auto stop = std::chrono::steady_clock::now();
  rec.cpu_ms = std::chrono::duration< double, std::milli >( stop - start )
    .count();
  counters.cpu_ms += rec.cpu_ms;

  if ( !ok && config_.verbose > 0 ) {
    G4cout << "[g4marley] MARLEY failed for Z=" << frag.GetZ_asInt() << " A="
      << frag.GetA_asInt() << " Ex=" << frag.GetExcitationEnergy()
      << " MeV after " << config_.max_attempts << " attempts (" << last_error
      << "); using native Geant4" << G4endl;
  }
  return ok;
}

int G4MarleyDeexEngine::choose_spin_parity( const G4Fragment& frag,
  int& twoJ, int& parity )
{
  const int pdg = marley_utils::get_nucleus_pid( frag.GetZ_asInt(),
    frag.GetA_asInt() );
  switch ( config_.spin_model ) {
    case G4MarleySpinModel::MarleyDefault:
      twoJ = -1;
      parity = 0;
      return 0;
    case G4MarleySpinModel::Fixed:
      twoJ = config_.fixed_twoJ;
      parity = config_.fixed_parity;
      return 3;
    case G4MarleySpinModel::Auto:
      if ( this->entrance_channel_spin_parity(frag, twoJ, parity) ) return 2;
      [[fallthrough]];
    case G4MarleySpinModel::LevelDensity:
    default:
      this->level_density_spin_parity( pdg, frag.GetExcitationEnergy(), twoJ,
        parity );
      return 1;
  }
}

void G4MarleyDeexEngine::level_density_spin_parity( int pdg, double Ex,
  int& twoJ, int& parity )
{
  const int A = marley_utils::get_particle_A( pdg );
  const int start_twoJ = ( A % 2 == 1 ) ? 1 : 0;

  // Weight each (J, Pi) by MARLEY's own level density at Ex
  auto& ld = gen_->get_structure_db().get_level_density_model( pdg );
  std::vector< double > weights;
  std::vector< std::pair<int, int> > spin_parities;
  for ( int tj = start_twoJ; tj <= config_.max_twoJ; tj += 2 ) {
    for ( int par : { 1, -1 } ) {
      double rho = ld.level_density( Ex, tj, marley::Parity(par) );
      if ( !std::isfinite(rho) || rho < 0. ) rho = 0.;
      weights.push_back( rho );
      spin_parities.emplace_back( tj, par );
    }
  }

  double sum = 0.;
  for ( double w : weights ) sum += w;
  if ( sum <= 0. ) {
    twoJ = start_twoJ;
    parity = 1;
    return;
  }

  std::discrete_distribution< size_t > dist( weights.begin(), weights.end() );
  size_t idx = gen_->sample_from_distribution( dist );
  twoJ = spin_parities.at( idx ).first;
  parity = spin_parities.at( idx ).second;
}

bool G4MarleyDeexEngine::entrance_channel_spin_parity( const G4Fragment& frag,
  int& twoJ, int& parity )
{
  // Find the interaction Geant4 is processing right now. This works without
  // any Geant4 modification because the de-excitation is called from inside
  // G4HadronicProcess::PostStepDoIt (or AtRestDoIt) of the current track.
  auto* evm = G4EventManager::GetEventManager();
  if ( !evm ) return false; // e.g. standalone use of G4ExcitationHandler
  auto* trm = evm->GetTrackingManager();
  if ( !trm ) return false;
  const G4Track* trk = trm->GetTrack();
  auto* stm = trm->GetSteppingManager();
  if ( !trk || !stm ) return false;
  const auto* hproc = dynamic_cast< const G4HadronicProcess* >(
    stm->GetfCurrentProcess() );
  if ( !hproc ) return false;
  const G4Nucleus* target = hproc->GetTargetNucleus();
  if ( !target ) return false;
  const int Zt = target->GetZ_asInt();
  const int At = target->GetA_asInt();

  // Projectile charge, mass number and 2*spin
  const int proj_pdg = trk->GetDefinition()->GetPDGEncoding();
  int zp = 0, ap = 0, two_s = 0;
  switch ( proj_pdg ) {
    case marley_utils::NEUTRON:  zp = 0; ap = 1; two_s = 1; break;
    case marley_utils::PROTON:   zp = 1; ap = 1; two_s = 1; break;
    case marley_utils::DEUTERON: zp = 1; ap = 2; two_s = 2; break;
    case marley_utils::TRITON:   zp = 1; ap = 3; two_s = 1; break;
    case marley_utils::HELION:   zp = 2; ap = 3; two_s = 1; break;
    case marley_utils::ALPHA:    zp = 2; ap = 4; two_s = 0; break;
    default: return false;
  }

  // Only a compound nucleus that has emitted nothing yet "remembers" the
  // entrance channel
  if ( frag.GetZ_asInt() != Zt + zp || frag.GetA_asInt() != At + ap ) {
    return false;
  }

  // Exact CM kinetic energy for a target at rest
  const double m_p = trk->GetDefinition()->GetPDGMass();
  const double M_t = G4NucleiProperties::GetNuclearMass( At, Zt );
  const double T_lab = trk->GetKineticEnergy();
  const double s = m_p*m_p + M_t*M_t + 2.*M_t*( T_lab + m_p );
  const double KE_cm = std::sqrt( s ) - m_p - M_t;
  if ( !(KE_cm > 0.) ) return false;

  entrance_.proj_pdg = proj_pdg;
  entrance_.proj_KE = T_lab;
  entrance_.target_pdg = marley_utils::get_nucleus_pid( Zt, At );

  int twoI = 0;
  marley::Parity Pt;
  marley::StructureDatabase::get_gs_spin_parity( Zt, At, twoI, Pt );
  const int pi_t = Pt ? 1 : -1;

  auto& om = gen_->get_structure_db().get_optical_model(
    marley_utils::get_nucleus_pid(Zt, At) );

  // sigma(J, Pi) ~ (2J + 1) sum_{l, j} T_lj, with j = l + s and J = j + I
  // (the common 1/((2I+1)(2s+1)) factor drops out after normalization)
  std::map< std::pair<int, int>, double > weights;
  constexpr int L_MAX = 15;
  for ( int l = 0; l <= L_MAX; ++l ) {
    const int pi = pi_t * ( l % 2 == 0 ? 1 : -1 );
    double sum_l = 0.;
    for ( int two_j = std::abs(2*l - two_s); two_j <= 2*l + two_s;
      two_j += 2 )
    {
      double T = om.transmission_coefficient( KE_cm, proj_pdg, two_j, l,
        two_s );
      if ( !std::isfinite(T) || T <= 0. ) continue;
      sum_l += T;
      for ( int tJ = std::abs(two_j - twoI); tJ <= two_j + twoI; tJ += 2 ) {
        weights[ {tJ, pi} ] += ( tJ + 1 ) * T;
      }
    }
    // Higher partial waves are negligible once T_l has died away
    if ( l > 2 && sum_l < 1e-8 ) break;
  }
  if ( weights.empty() ) return false;

  std::vector< double > w;
  std::vector< std::pair<int, int> > sp;
  for ( const auto& kv : weights ) {
    sp.push_back( kv.first );
    w.push_back( kv.second );
  }
  std::discrete_distribution< size_t > dist( w.begin(), w.end() );
  size_t idx = gen_->sample_from_distribution( dist );
  twoJ = sp.at( idx ).first;
  parity = sp.at( idx ).second;
  return true;
}

bool G4MarleyDeexEngine::try_decay_once( const G4Fragment& frag,
  G4MarleyDeexRecord& rec )
{
  const int Z = frag.GetZ_asInt();
  const int A = frag.GetA_asInt();
  const double Ex = frag.GetExcitationEnergy();
  const int pdg = marley_utils::get_nucleus_pid( Z, A );

  // MARLEY works with atomic masses. A Geant4 fragment is a bare nucleus,
  // so give it ionic charge Z (MARLEY then subtracts Z electron masses).
  const auto& mt = marley::MassTable::Instance();
  const double me = mt.get_particle_mass( marley_utils::ELECTRON );
  const double M_gs_marley = mt.get_atomic_mass( pdg ) - Z*me;
  const double M_marley = M_gs_marley + Ex;

  int twoJ = -1, parity = 0;
  rec.spin_source = this->choose_spin_parity( frag, twoJ, parity );
  rec.twoJ = twoJ;
  rec.parity = parity;
  if ( rec.spin_source == 2 ) {
    rec.proj_pdg = entrance_.proj_pdg;
    rec.proj_KE = entrance_.proj_KE;
    rec.target_pdg = entrance_.target_pdg;
  }

  // Build a one-particle event: the excited nucleus at rest, flagged as an
  // undecayed residue (NuHepMC status 27) so NucleusDecayer will pick it up.
  // There is no MARLEY GenRunInfo attached, so NucleusDecayer takes its
  // "externally produced event" path (Ex inferred from the mass, continuum
  // Hauser-Feshbach start).
  HepMC3::GenEvent ev( HepMC3::Units::MEV, HepMC3::Units::CM );
  auto prod_vtx = std::make_shared< HepMC3::GenVertex >();
  auto residue = marley_hepmc3::make_particle( pdg, 0., 0., 0., M_marley,
    marley_hepmc3::NUHEPMC_UNDECAYED_RESIDUE_STATUS, M_marley );
  prod_vtx->add_particle_out( residue );
  ev.add_vertex( prod_vtx );

  // Attributes can only be attached once the particle belongs to an event
  marley_hepmc3::set_particle_charge( *residue, Z );
  if ( twoJ >= 0 ) {
    // Honoured by the g4-deexcitation branch of NucleusDecayer
    residue->add_attribute( "twoJ",
      std::make_shared< HepMC3::IntAttribute >(twoJ) );
    residue->add_attribute( "parity",
      std::make_shared< HepMC3::IntAttribute >(parity) );
  }

  marley::NucleusDecayer decayer;
  decayer.process_event( ev, *gen_ );

  // ---- Read back the cascade (for bookkeeping/annotation only) ----
  for ( const auto& vtx : ev.vertices() ) {
    const int st = vtx->status();
    if ( st != marley_hepmc3::NUHEPMC_HF_DECAY_VERTEX
      && st != marley_hepmc3::NUHEPMC_GAMMA_DECAY_VERTEX ) continue;
    if ( vtx->particles_in().empty() || vtx->particles_out().size() < 2 ) {
      continue;
    }
    const auto& parent = vtx->particles_in().front();
    const auto& emitted = vtx->particles_out().at( 0 );
    const auto& daughter = vtx->particles_out().at( 1 );

    G4MarleyDeexStep step;
    step.type = st;
    step.parent_pdg = parent->pid();
    if ( parent == residue ) {
      step.parent_Ex = Ex;
      step.parent_twoJ = twoJ;
      step.parent_parity = parity;
    }
    else {
      step.parent_Ex = double_attr( parent, "Ex", -1. );
      step.parent_twoJ = int_attr( parent, "twoJ", -1 );
      step.parent_parity = int_attr( parent, "parity", 0 );
    }
    step.emitted_pdg = emitted->pid();
    step.emitted_KE = emitted->momentum().e() - emitted->generated_mass();
    if ( emitted->pid() == marley_utils::PHOTON ) {
      step.emitted_KE = emitted->momentum().e();
    }
    step.daughter_pdg = daughter->pid();
    step.daughter_Ex = double_attr( daughter, "Ex", -1. );
    step.daughter_twoJ = int_attr( daughter, "twoJ", -1 );
    step.daughter_parity = int_attr( daughter, "parity", 0 );
    step.time_ns = vtx->position().t() / C_CM_PER_NS;
    rec.steps.push_back( step );
  }

  // ---- Collect the final state: every particle without an end vertex ----
  struct Item {
    int pdg;
    G4ThreeVector p;
    double m; // Geant4 mass
    double t_ns;
    bool is_residual;
    int emitter_pdg; // nucleus at the decay vertex that produced it
  };
  std::vector< Item > items;
  int residual_index = -1;

  for ( const auto& p : ev.particles() ) {
    if ( p->end_vertex() ) continue;
    if ( p == residue ) return false; // nothing happened
    const int fpdg = p->pid();
    const double m = g4_mass( fpdg );
    if ( m < 0. ) {
      throw std::runtime_error( "unexpected MARLEY product PDG "
        + std::to_string(fpdg) );
    }
    const auto& mom = p->momentum();
    Item it;
    it.pdg = fpdg;
    it.p = G4ThreeVector( mom.px(), mom.py(), mom.pz() );
    it.m = m;
    it.t_ns = p->production_vertex()
      ? p->production_vertex()->position().t() / C_CM_PER_NS : 0.;
    it.is_residual = ( p->status()
      == marley_hepmc3::NUHEPMC_INTERMEDIATE_RESIDUE_STATUS );
    it.emitter_pdg = 0;
    if ( !it.is_residual && p->production_vertex()
      && !p->production_vertex()->particles_in().empty() )
    {
      it.emitter_pdg = p->production_vertex()->particles_in().front()->pid();
    }
    if ( it.is_residual ) {
      // Excitation MARLEY left in the final nucleus (e.g., the 1 keV
      // continuum cut-off). It is absorbed by the rescaling below.
      const int q = marley_hepmc3::get_particle_charge( *p );
      rec.ex_left = p->generated_mass() - ( mt.get_atomic_mass(fpdg) - q*me );
      residual_index = static_cast< int >( items.size() );
    }
    items.push_back( it );
  }

  if ( items.size() < 2 ) return false;

  // Fall back to "heaviest product" if no residual status was found
  if ( residual_index < 0 ) {
    int best_A = -1;
    for ( size_t i = 0; i < items.size(); ++i ) {
      if ( !marley_utils::is_ion(items[i].pdg) ) continue;
      int a = marley_utils::get_particle_A( items[i].pdg );
      if ( a > best_A ) { best_A = a; residual_index = static_cast<int>( i ); }
    }
    if ( residual_index < 0 ) return false;
  }
  items[ residual_index ].emitter_pdg = 0;

  // ---- Energy-conserving rescale with Geant4 masses ----
  // Solve sum_i sqrt( (lambda*p_i)^2 + m_i^2 ) = M_G4 for lambda by Newton's
  // method. The left side increases monotonically with lambda.
  const double M_G4 = frag.GetMomentum().m();
  double sum_m = 0.;
  for ( const auto& it : items ) sum_m += it.m;
  if ( sum_m >= M_G4 ) return false; // not enough energy with G4 masses

  double lambda = 1.;
  for ( int iter = 0; iter < 100; ++iter ) {
    double f = -M_G4, df = 0.;
    for ( const auto& it : items ) {
      const double p2 = it.p.mag2();
      const double E = std::sqrt( lambda*lambda*p2 + it.m*it.m );
      f += E;
      if ( E > 0. ) df += lambda * p2 / E;
    }
    if ( std::abs(f) < E_TOLERANCE ) break;
    if ( df <= 0. ) return false;
    double next = lambda - f / df;
    // Keep the iteration in the physical region
    lambda = ( next > 0. ) ? next : 0.5*lambda;
  }
  rec.lambda = lambda;

  // ---- Boost to the frame Geant4 gave us and store ----
  const G4ThreeVector beta = frag.GetMomentum().boostVector();
  auto make_product = [&]( const Item& it ) {
    G4ThreeVector p = lambda * it.p;
    G4LorentzVector lv( p, std::sqrt(p.mag2() + it.m*it.m) );
    lv.boost( beta );
    G4MarleyDeexProduct prod;
    prod.pdg = it.pdg;
    prod.p4 = lv;
    prod.time_ns = it.t_ns;
    prod.emitter_pdg = it.emitter_pdg;
    return prod;
  };

  for ( size_t i = 0; i < items.size(); ++i ) {
    if ( static_cast<int>(i) == residual_index ) continue;
    rec.products.push_back( make_product(items[i]) );
  }
  rec.products.push_back( make_product(items[residual_index]) );

  return true;
}
