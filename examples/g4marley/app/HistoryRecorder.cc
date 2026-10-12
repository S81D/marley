/// @file
/// @brief Writes the full Geant4 particle history of each event to HepMC3

#include <cmath>
#include <iostream>
#include <sstream>

#include "G4Event.hh"
#include "G4HadronicProcess.hh"
#include "G4Nucleus.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicalConstants.hh"
#include "G4Run.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4Version.hh"

#include "HepMC3/Attribute.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"

#include "AppActions.hh"
#include "HistoryRecorder.hh"

namespace {

  // c in mm/ns: HepMC3 stores time as c*t in the length unit (mm here)
  constexpr double C_MM_PER_NS = 299.792458;

  // Vertex status codes (see HistoryRecorder.hh)
  int vertex_status( G4ProcessType type ) {
    switch ( type ) {
      case fHadronic: return 11;
      case fDecay: return 12;
      case fElectromagnetic: return 13;
      default: return 14;
    }
  }

  std::string process_name( const G4VProcess* p ) {
    return p ? std::string( p->GetProcessName() ) : std::string( "none" );
  }

  bool is_em_lepton_or_photon( int pdg ) {
    return pdg == 22 || pdg == 11 || pdg == -11;
  }

}

HistoryRecorder::HistoryRecorder( const AppOptions& opts ) : opts_( opts ) {
  if ( !opts_.deex_log_file.empty() ) deex_log_.open( opts_.deex_log_file );
  if ( !opts_.had_csv_file.empty() ) {
    had_csv_.open( opts_.had_csv_file );
    had_csv_ << "evt,vtx,proc,proj_pdg,proj_KE,tgt_Z,tgt_A,sec_pdg,sec_KE,deex\n";
  }
}

HistoryRecorder::~HistoryRecorder() {
  if ( writer_ ) writer_->close();
}

HepMC3::GenParticlePtr HistoryRecorder::make_particle( const G4Track* trk,
  int status ) const
{
  const G4ThreeVector p = trk->GetMomentum();
  HepMC3::FourVector p4( p.x()/MeV, p.y()/MeV, p.z()/MeV,
    trk->GetTotalEnergy()/MeV );
  auto part = std::make_shared< HepMC3::GenParticle >( p4,
    trk->GetDefinition()->GetPDGEncoding(), status );
  part->set_generated_mass( trk->GetDefinition()->GetPDGMass()/MeV );
  return part;
}

bool HistoryRecorder::keep_secondary( const G4Track* sec ) const {
  if ( opts_.record_min_ekin_MeV <= 0. ) return true;
  int pdg = sec->GetDefinition()->GetPDGEncoding();
  if ( !is_em_lepton_or_photon(pdg) ) return true;
  return sec->GetKineticEnergy()/MeV >= opts_.record_min_ekin_MeV;
}

void HistoryRecorder::push_point( OpenTrack& o, const G4ThreeVector& pos,
  double t_ns )
{
  o.traj.push_back( pos.x()/mm );
  o.traj.push_back( pos.y()/mm );
  o.traj.push_back( pos.z()/mm );
  o.traj.push_back( t_ns * C_MM_PER_NS );
}

void HistoryRecorder::finish_segment( OpenTrack& o, bool has_end_vertex ) {
  // Segments ending in a vertex only need the attribute if they have
  // intermediate points (the end point is the vertex itself).
  if ( !o.particle ) return;
  bool worth_it = !o.traj.empty()
    && ( !has_end_vertex || o.traj.size() > 4 );
  if ( worth_it && ( opts_.store_trajectories || !has_end_vertex ) ) {
    o.particle->add_attribute( "traj",
      std::make_shared< HepMC3::VectorDoubleAttribute >(o.traj) );
  }
  o.traj.clear();
}

void HistoryRecorder::begin_event( const G4Event* event ) {
  if ( !run_info_ ) {
    run_info_ = std::make_shared< HepMC3::GenRunInfo >();
    run_info_->tools().push_back( { "Geant4", G4Version, "" } );
    run_info_->tools().push_back( { "MARLEY", "g4-deexcitation branch",
      "nuclear de-excitation" } );
    run_info_->tools().push_back( { "g4marley", "0.1",
      "MARLEY de-excitation inside Geant4 (test application)" } );
    run_info_->add_attribute( "conventions",
      std::make_shared< HepMC3::StringAttribute >(
      "units MeV,mm; vertex t = c*t (mm); one particle per track segment;"
      " particle attrs tid,ptid,cproc,eproc,traj(x,y,z,ct...),deex;"
      " vertex attrs proc,marley(JSON); primaries attached to root vertex"
      " (event position); vertex status 11 hadronic"
      " 12 decay 13 EM 14 other; particle status 4 primary 2 interacted"
      " 1 no end vertex" ) );
    for ( const auto& kv : run_attrs_ ) {
      run_info_->add_attribute( kv.first,
        std::make_shared< HepMC3::StringAttribute >(kv.second) );
    }
    if ( opts_.hepmc_file != "none" ) {
      writer_ = std::make_unique< HepMC3::WriterAscii >( opts_.hepmc_file,
        run_info_ );
    }
  }

  current_event_ = event->GetEventID();
  event_ = std::make_unique< HepMC3::GenEvent >( run_info_,
    HepMC3::Units::MEV, HepMC3::Units::MM );
  event_->set_event_number( current_event_ );
  open_.clear();
  evt_marley_accepted_ = 0;
  evt_marley_discarded_ = 0;
  // Anything left over from a previous event is stale
  G4MarleyDeexRecorder::Instance().clear_pending();
}

void HistoryRecorder::end_event( const G4Event* ) {
  const auto& gen = CurrentGenEventInfo();
  if ( !opts_.hepmc_in_file.empty() && gen.input_index >= 0 ) {
    event_->add_attribute( "gen_input_index",
      std::make_shared< HepMC3::IntAttribute >(static_cast<int>( gen.input_index )) );
    event_->add_attribute( "gen_input_event",
      std::make_shared< HepMC3::IntAttribute >(gen.input_event_number) );
    for ( const auto& kv : gen.attrs ) {
      event_->add_attribute( kv.first, std::make_shared< HepMC3::StringAttribute >(
        kv.second ) );
    }
  }
  event_->add_attribute( "marley_calls",
    std::make_shared< HepMC3::IntAttribute >(evt_marley_accepted_) );
  event_->add_attribute( "marley_discarded",
    std::make_shared< HepMC3::IntAttribute >(evt_marley_discarded_) );
  if ( writer_ ) writer_->write_event( *event_ );
  event_.reset();
  open_.clear();
}

void HistoryRecorder::pre_track( const G4Track* track ) {
  OpenTrack o;
  const int tid = track->GetTrackID();
  auto* info = dynamic_cast< HistoryTrackInfo* >(
    track->GetUserInformation() );

  if ( track->GetParentID() == 0 ) {
    // Primary: attach it to the event's root vertex (HepMC3 convention for
    // particles without a production vertex; the ASCII format cannot hold a
    // vertex with no incoming particles). The root vertex position is the
    // gun position.
    if ( event_->vertices().empty() && event_->particles().empty() ) {
      const G4ThreeVector pos = track->GetPosition();
      event_->shift_position_to( HepMC3::FourVector( pos.x()/mm, pos.y()/mm,
        pos.z()/mm, track->GetGlobalTime()/ns * C_MM_PER_NS ) );
    }
    o.particle = make_particle( track, 4 );
    event_->add_particle( o.particle );
    o.particle->add_attribute( "cproc",
      std::make_shared< HepMC3::StringAttribute >("primary") );
    o.recorded = true;
  }
  else if ( info ) {
    o.particle = info->particle;
    o.recorded = info->recorded && o.particle;
  }

  if ( o.recorded ) {
    o.particle->add_attribute( "tid",
      std::make_shared< HepMC3::IntAttribute >(tid) );
    o.particle->add_attribute( "ptid",
      std::make_shared< HepMC3::IntAttribute >(track->GetParentID()) );
  }
  o.last_pos = track->GetPosition();
  o.last_dir = track->GetMomentumDirection();
  open_[ tid ] = std::move( o );
}

void HistoryRecorder::post_track( const G4Track* track ) {
  auto it = open_.find( track->GetTrackID() );
  if ( it == open_.end() ) return;
  OpenTrack& o = it->second;
  if ( o.recorded && o.particle && !o.particle->end_vertex() ) {
    push_point( o, track->GetPosition(), track->GetGlobalTime()/ns );
    finish_segment( o, false );
    // Process that ended the track (Transportation = left the world)
    const G4Step* last = track->GetStep();
    const G4VProcess* proc = last ? last->GetPostStepPoint()
      ->GetProcessDefinedStep() : nullptr;
    o.particle->add_attribute( "eproc",
      std::make_shared< HepMC3::StringAttribute >(process_name(proc)) );
  }
  open_.erase( it );
}

std::vector< const G4MarleyDeexRecord* > HistoryRecorder::match_marley(
  const std::vector< G4MarleyDeexRecord >& recs,
  const std::vector< const G4Track* >& secs,
  std::vector< long >& matched_ids ) const
{
  std::vector< const G4MarleyDeexRecord* > accepted;
  matched_ids.assign( secs.size(), 0 );

  // Latest record first: when Bertini retries a cascade, only the last
  // attempt's products survive into the secondaries.
  for ( auto r = recs.rbegin(); r != recs.rend(); ++r ) {
    std::vector< long > trial = matched_ids;
    bool all_ok = true;
    for ( size_t ip = 0; ip < r->products.size(); ++ip ) {
      const auto& prod = r->products[ ip ];
      const bool is_residual = ( ip + 1 == r->products.size() );
      double ke = prod.p4.e() - prod.p4.m();
      if ( prod.pdg == 22 ) ke = prod.p4.e();
      // Kinetic energies are invariant under the rotation Bertini applies
      // between its internal frame and the lab
      double tol = 1e-5*ke + 1e-7;
      int best = -1;
      double best_diff = tol;
      for ( size_t is = 0; is < secs.size(); ++is ) {
        if ( trial[is] != 0 ) continue;
        if ( secs[is]->GetDefinition()->GetPDGEncoding() != prod.pdg ) continue;
        double diff = std::abs( secs[is]->GetKineticEnergy()/MeV - ke );
        if ( diff <= best_diff ) { best_diff = diff; best = static_cast<int>(is); }
      }
      if ( best >= 0 ) trial[ best ] = r->call_id;
      else if ( !is_residual ) { all_ok = false; break; }
    }
    if ( all_ok ) {
      matched_ids = trial;
      accepted.push_back( &(*r) );
    }
  }
  return accepted;
}

void HistoryRecorder::write_deex_log( const G4Step* step,
  const std::vector< G4MarleyDeexRecord >& recs,
  const std::vector< const G4MarleyDeexRecord* >& accepted )
{
  if ( !deex_log_.is_open() ) return;
  const G4Track* trk = step->GetTrack();
  for ( const auto& r : recs ) {
    bool acc = false;
    for ( const auto* a : accepted ) if ( a == &r ) acc = true;
    deex_log_ << "{\"evt\":" << current_event_ << ",\"tid\":"
      << trk->GetTrackID() << ",\"proj\":"
      << trk->GetDefinition()->GetPDGEncoding() << ",\"proj_KE\":"
      << step->GetPreStepPoint()->GetKineticEnergy()/MeV << ",\"proc\":\""
      << process_name( step->GetPostStepPoint()->GetProcessDefinedStep() )
      << "\",\"accepted\":" << ( acc ? "true" : "false" )
      << ",\"n_calls_step\":" << recs.size()
      << ",\"rec\":" << r.to_json() << "}\n";
  }
}

void HistoryRecorder::step( const G4Step* step ) {
  const G4Track* trk = step->GetTrack();
  auto it = open_.find( trk->GetTrackID() );
  if ( it == open_.end() ) return;
  OpenTrack& o = it->second;

  const auto* sec_vec = step->GetSecondaryInCurrentStep();
  std::vector< const G4Track* > secs( sec_vec->begin(), sec_vec->end() );

  const G4StepPoint* post = step->GetPostStepPoint();
  const G4VProcess* proc = post->GetProcessDefinedStep();
  const std::string pname = process_name( proc );
  const G4ProcessType ptype = proc ? proc->GetProcessType() : fNotDefined;

  // ---- MARLEY bookkeeping (independent of what is written) ----
  auto recs = G4MarleyDeexRecorder::Instance().take_pending();
  std::vector< long > matched_ids;
  auto accepted = match_marley( recs, secs, matched_ids );
  if ( ptype == fHadronic ) {
    ++calls_per_step_[ static_cast<int>( recs.size() ) ];
    if ( secs.empty() && trk->GetTrackStatus() == fAlive
      && std::abs( post->GetKineticEnergy()
      - step->GetPreStepPoint()->GetKineticEnergy() ) < 1e-12
      && post->GetMomentumDirection() == step->GetPreStepPoint()
      ->GetMomentumDirection() )
    {
      ++run_no_interaction_steps_;
    }
  }
  evt_marley_accepted_ += accepted.size();
  evt_marley_discarded_ += recs.size() - accepted.size();
  run_marley_accepted_ += accepted.size();
  run_marley_discarded_ += recs.size() - accepted.size();
  write_deex_log( step, recs, accepted );

  // ---- Compact record of non-elastic hadronic interactions ----
  const bool is_elastic = pname.find( "Elastic" ) != std::string::npos
    && pname.find( "Inelastic" ) == std::string::npos;
  if ( had_csv_.is_open() && ptype == fHadronic && !is_elastic ) {
    ++had_vertex_counter_;
    int tZ = 0, tA = 0;
    if ( const auto* hp = dynamic_cast< const G4HadronicProcess* >( proc ) ) {
      if ( const G4Nucleus* nuc = hp->GetTargetNucleus() ) {
        tZ = nuc->GetZ_asInt();
        tA = nuc->GetA_asInt();
      }
    }
    std::ostringstream head;
    head << current_event_ << ',' << had_vertex_counter_ << ',' << pname << ','
      << trk->GetDefinition()->GetPDGEncoding() << ','
      << step->GetPreStepPoint()->GetKineticEnergy()/MeV << ',' << tZ << ','
      << tA << ',';
    if ( secs.empty() ) had_csv_ << head.str() << "0,0,0\n";
    for ( size_t i = 0; i < secs.size(); ++i ) {
      had_csv_ << head.str() << secs[i]->GetDefinition()->GetPDGEncoding()
        << ',' << secs[i]->GetKineticEnergy()/MeV << ',' << matched_ids[i]
        << '\n';
    }
  }

  // ---- Unrecorded subtree: propagate the flag and stop ----
  if ( !o.recorded ) {
    for ( const auto* s : secs ) {
      const_cast< G4Track* >( s )->SetUserInformation(
        new HistoryTrackInfo(nullptr, false) );
    }
    return;
  }

  const G4TrackStatus status = trk->GetTrackStatus();
  const bool alive = ( status == fAlive || status == fStopButAlive
    || status == fSuspend );
  const bool charged = ( trk->GetDefinition()->GetPDGCharge() != 0. );
  const bool transport = ( ptype == fTransportation || ptype == fNotDefined );

  std::vector< bool > keep( secs.size() );
  int n_keep = 0;
  for ( size_t i = 0; i < secs.size(); ++i ) {
    keep[i] = keep_secondary( secs[i] );
    if ( keep[i] ) ++n_keep;
  }

  // An interaction vertex: something was produced, or a hadronic/decay
  // process acted, or a neutral particle had any non-transport interaction
  // (neutral tracks are straight lines between such points).
  // A vertex also needs at least one outgoing particle: the HepMC3 ASCII
  // format cannot store a vertex without outgoing particles. A track that
  // just ends (absorbed with nothing recorded, stopped, ...) is handled in
  // post_track() instead, with the ending process in the "eproc" attribute.
  const bool make_vertex = ( ( n_keep > 0 ) || ptype == fHadronic
    || ptype == fDecay || ( !charged && !transport ) )
    && ( n_keep > 0 || alive );

  if ( !make_vertex ) {
    // Not a vertex: maybe store a trajectory point (charged tracks only)
    if ( opts_.store_trajectories && charged && alive ) {
      const G4ThreeVector pos = post->GetPosition();
      const G4ThreeVector dir = post->GetMomentumDirection();
      double cosang = dir.dot( o.last_dir );
      if ( cosang < std::cos(5.*deg) || (pos - o.last_pos).mag() > 20.*mm ) {
        if ( o.traj.size() < 4*400 ) {
          push_point( o, pos, post->GetGlobalTime()/ns );
          o.last_pos = pos;
          o.last_dir = dir;
        }
      }
    }
    for ( const auto* s : secs ) {
      const_cast< G4Track* >( s )->SetUserInformation(
        new HistoryTrackInfo(nullptr, false) );
    }
    return;
  }

  // ---- Create the vertex ----
  const G4ThreeVector pos = post->GetPosition();
  const double t_ns = post->GetGlobalTime()/ns;
  auto vtx = std::make_shared< HepMC3::GenVertex >( HepMC3::FourVector(
    pos.x()/mm, pos.y()/mm, pos.z()/mm, t_ns * C_MM_PER_NS ) );
  vtx->set_status( vertex_status(ptype) );
  event_->add_vertex( vtx );
  vtx->add_particle_in( o.particle );
  if ( o.particle->status() != 4 ) o.particle->set_status( 2 );
  finish_segment( o, true );
  vtx->add_attribute( "proc", std::make_shared< HepMC3::StringAttribute >(
    pname ) );

  for ( size_t i = 0; i < secs.size(); ++i ) {
    if ( !keep[i] ) {
      const_cast< G4Track* >( secs[i] )->SetUserInformation(
        new HistoryTrackInfo(nullptr, false) );
      continue;
    }
    auto p = make_particle( secs[i], 1 );
    vtx->add_particle_out( p );
    p->add_attribute( "cproc", std::make_shared< HepMC3::StringAttribute >(
      process_name(secs[i]->GetCreatorProcess()) ) );
    if ( matched_ids[i] != 0 ) {
      p->add_attribute( "deex", std::make_shared< HepMC3::IntAttribute >(
        static_cast<int>( matched_ids[i] ) ) );
    }
    const_cast< G4Track* >( secs[i] )->SetUserInformation(
      new HistoryTrackInfo(p, true) );
  }

  // The track survives: continue it as a new segment from this vertex
  if ( alive ) {
    auto p = make_particle( trk, 1 );
    vtx->add_particle_out( p );
    p->add_attribute( "tid", std::make_shared< HepMC3::IntAttribute >(
      trk->GetTrackID() ) );
    p->add_attribute( "ptid", std::make_shared< HepMC3::IntAttribute >(
      trk->GetParentID() ) );
    p->add_attribute( "cproc", std::make_shared< HepMC3::StringAttribute >(
      "continued" ) );
    o.particle = p;
    o.last_pos = pos;
    o.last_dir = post->GetMomentumDirection();
  }

  if ( !accepted.empty() ) {
    std::string js = "[";
    for ( size_t i = 0; i < accepted.size(); ++i ) {
      if ( i ) js += ",";
      js += accepted[i]->to_json();
    }
    js += "]";
    vtx->add_attribute( "marley", std::make_shared< HepMC3::StringAttribute >(
      js ) );
  }
}

void HistoryRecorder::end_of_run( const G4Run* run ) {
  std::cout << "\n==== g4marley run summary (" << run->GetNumberOfEvent()
    << " events) ====\n"
    << "  MARLEY calls whose products became Geant4 secondaries: "
    << run_marley_accepted_ << '\n'
    << "  MARLEY calls discarded (e.g. Bertini retries):         "
    << run_marley_discarded_ << '\n'
    << "  hadronic steps with no effect (NoInteraction):          "
    << run_no_interaction_steps_ << '\n'
    << "  MARLEY calls per hadronic step (calls: steps):";
  for ( const auto& kv : calls_per_step_ ) {
    std::cout << "  " << kv.first << ": " << kv.second;
  }
  std::cout << '\n';
  G4MarleyDeexRecorder::Instance().print_summary( std::cout );
  if ( writer_ ) {
    writer_->close();
    writer_.reset();
  }
  if ( deex_log_.is_open() ) deex_log_.close();
  if ( had_csv_.is_open() ) had_csv_.close();
}
