/// @file
/// @brief G4VEvaporation implementation that hands nuclei to MARLEY
//
// This file is part of the MARLEY Geant4 de-excitation interface (g4marley).
// It is distributed under the same terms as MARLEY (GPL v3).

#include "G4Evaporation.hh"
#include "G4Fragment.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4ios.hh"

#include "marley/marley_utils.hh"

#include "G4MarleyDeexEngine.hh"
#include "G4MarleyDeexRecorder.hh"
#include "G4MarleyEvaporation.hh"

G4MarleyEvaporation::G4MarleyEvaporation(
  std::shared_ptr< G4MarleyDeexEngine > engine )
  : G4VEvaporation(), engine_( std::move(engine) ),
  fallback_( new G4Evaporation() )
{
  // G4Evaporation() creates its own G4PhotonEvaporation. Mirror the pointer
  // so G4ExcitationHandler::SetEvaporation() adopts that same object.
  thePhotonEvaporation = fallback_->GetPhotonEvaporation();
}

G4MarleyEvaporation::~G4MarleyEvaporation() {
  // The fallback owns the photon evaporation; stop ~G4VEvaporation() from
  // deleting it a second time.
  thePhotonEvaporation = nullptr;
  delete fallback_;
}

void G4MarleyEvaporation::SetPhotonEvaporation( G4VEvaporationChannel* ptr ) {
  fallback_->SetPhotonEvaporation( ptr );
  thePhotonEvaporation = fallback_->GetPhotonEvaporation();
}

void G4MarleyEvaporation::sync_fallback() {
  // SetFermiBreakUp() is non-virtual, so G4ExcitationHandler only sets it on
  // this object. Forward it before the fallback can use it.
  if ( theFBU && fallback_->GetFermiBreakUp() != theFBU ) {
    fallback_->SetFermiBreakUp( theFBU );
  }
}

void G4MarleyEvaporation::InitialiseChannels() {
  sync_fallback();
  fallback_->InitialiseChannels();
}

void G4MarleyEvaporation::BreakFragment( G4FragmentVector* products,
  G4Fragment* nucleus )
{
  sync_fallback();

  auto& recorder = G4MarleyDeexRecorder::Instance();
  auto& counters = recorder.counters();

  auto skip = engine_->check( *nucleus );
  if ( skip != G4MarleyDeexEngine::Skip::None ) {
    switch ( skip ) {
      case G4MarleyDeexEngine::Skip::LowEx: ++counters.native_low_Ex; break;
      case G4MarleyDeexEngine::Skip::HighEx: ++counters.native_high_Ex; break;
      case G4MarleyDeexEngine::Skip::Light: ++counters.native_light; break;
      case G4MarleyDeexEngine::Skip::Disabled: ++counters.native_disabled;
        break;
      default: break;
    }
    fallback_->BreakFragment( products, nucleus );
    return;
  }

  G4MarleyDeexRecord rec;
  if ( !engine_->decay(*nucleus, rec) ) {
    ++counters.marley_failed;
    fallback_->BreakFragment( products, nucleus );
    return;
  }
  ++counters.marley_calls;
  rec.call_id = recorder.next_call_id();

  const double t0 = nucleus->GetCreationTime();

  // All products except the last one (the residual) become new fragments
  for ( size_t i = 0; i + 1 < rec.products.size(); ++i ) {
    const auto& prod = rec.products[i];
    G4Fragment* frag = nullptr;
    if ( prod.pdg == marley_utils::PHOTON ) {
      frag = new G4Fragment( prod.p4, G4Gamma::Gamma() );
    }
    else if ( prod.pdg == marley_utils::ELECTRON ) {
      frag = new G4Fragment( prod.p4, G4Electron::Electron() );
    }
    else if ( prod.pdg == marley_utils::NEUTRON ) {
      frag = new G4Fragment( 1, 0, prod.p4 );
    }
    else if ( prod.pdg == marley_utils::PROTON ) {
      frag = new G4Fragment( 1, 1, prod.p4 );
    }
    else {
      frag = new G4Fragment( marley_utils::get_particle_A(prod.pdg),
        marley_utils::get_particle_Z(prod.pdg), prod.p4 );
    }
    frag->SetCreationTime( t0 + prod.time_ns );
    products->push_back( frag );
  }

  // The input fragment becomes the (ground-state) residual nucleus
  const auto& res = rec.products.back();
  nucleus->SetZandA_asInt( marley_utils::get_particle_Z(res.pdg),
    marley_utils::get_particle_A(res.pdg) );
  nucleus->SetNumberOfHoles( 0, 0 );
  nucleus->SetNumberOfExcitedParticle( 0, 0 );
  nucleus->SetMomentum( res.p4 ); // recomputes Ex (= 0 up to round-off)
  nucleus->SetCreationTime( t0 + res.time_ns );
  nucleus->SetSpin( 0. );

  if ( engine_->config().verbose > 0 ) {
    G4cout << "[g4marley] #" << rec.call_id << " " << rec.chain_string()
      << "  (" << rec.cpu_ms << " ms)" << G4endl;
  }

  if ( engine_->config().record ) recorder.add( std::move(rec) );
}
