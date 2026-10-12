/// @file
/// @brief Geometry, primary generator and user actions for g4marley

#include <cmath>

#include <cstdlib>

#include "G4Box.hh"
#include "G4Event.hh"
#include "G4IonTable.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PVPlacement.hh"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/ReaderAscii.h"

#include "AppActions.hh"
#include "HistoryRecorder.hh"

G4VPhysicalVolume* DetectorConstruction::Construct() {
  auto* lar = G4NistManager::Instance()->FindOrBuildMaterial( "G4_lAr" );
  auto* box = new G4Box( "World", half_m_*m, half_m_*m, half_m_*m );
  auto* lv = new G4LogicalVolume( box, lar, "World" );
  return new G4PVPlacement( nullptr, G4ThreeVector(), lv, "World", nullptr,
    false, 0 );
}

GenEventInfo& CurrentGenEventInfo() {
  static thread_local GenEventInfo info;
  return info;
}

PrimaryGeneratorAction::PrimaryGeneratorAction( const AppOptions& opts )
  : gun_( new G4ParticleGun(1) ), isotropic_( opts.gun_isotropic ),
  species_( opts.hepmc_in_species )
{
  if ( !opts.hepmc_in_file.empty() ) {
    reader_ = std::make_unique< HepMC3::ReaderAscii >( opts.hepmc_in_file );
    if ( reader_->failed() ) {
      G4cerr << "Cannot read HepMC3 input " << opts.hepmc_in_file << G4endl;
      std::exit( 3 );
    }
  }
  auto* def = G4ParticleTable::GetParticleTable()->FindParticle(
    opts.gun_particle );
  if ( def ) gun_->SetParticleDefinition( def );
  gun_->SetParticleEnergy( opts.gun_energy_MeV * MeV );
  gun_->SetParticlePosition( G4ThreeVector() );
  gun_->SetParticleMomentumDirection( G4ThreeVector(0., 0., 1.) );
  gun_->SetParticleTime( 0. );
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() { delete gun_; }

bool PrimaryGeneratorAction::keep_species( int pdg ) const {
  const int a = std::abs( pdg );
  if ( species_ == "all" ) return true;
  if ( species_ == "nucleons" ) return a == 2112 || a == 2212;
  return a < 11 || a > 16;                     // no-leptons
}

const G4ParticleDefinition* PrimaryGeneratorAction::definition( int pdg ) const {
  if ( pdg > 1000000000 ) {
    // Nucleus (10LZZZAAAI); injected in its ground state
    const int Z = ( pdg / 10000 ) % 1000;
    const int A = ( pdg / 10 ) % 1000;
    if ( Z < 1 || A < 1 || A < Z ) return nullptr;
    return G4IonTable::GetIonTable()->GetIon( Z, A, 0. );
  }
  // GENIE pseudo-particles (rootino, bindino, ...) are not found and skipped
  return G4ParticleTable::GetParticleTable()->FindParticle( pdg );
}

bool PrimaryGeneratorAction::generate_from_hepmc( G4Event* event ) {
  HepMC3::GenEvent in;
  // Skip input events with nothing to inject (e.g. no nucleons)
  while ( true ) {
    reader_->read_event( in );
    if ( reader_->failed() ) return false;
    in.set_units( HepMC3::Units::MEV, HepMC3::Units::MM );
    const long index = n_read_++;
    auto* vtx = new G4PrimaryVertex( gun_->GetParticlePosition(), 0. );
    int n = 0;
    for ( const auto& p : in.particles() ) {
      if ( p->status() != 1 || !keep_species(p->pid()) ) continue;
      const G4ParticleDefinition* def = definition( p->pid() );
      if ( !def ) continue;
      const auto& m = p->momentum();
      vtx->SetPrimary( new G4PrimaryParticle(def, m.px()*MeV, m.py()*MeV,
        m.pz()*MeV) );
      ++n;
    }
    if ( n == 0 ) { delete vtx; continue; }
    event->AddPrimaryVertex( vtx );
    auto& info = CurrentGenEventInfo();
    info.input_index = index;
    info.input_event_number = in.event_number();
    info.attrs.clear();
    for ( const auto& name : in.attribute_names() ) {
      info.attrs[ name ] = in.attribute_as_string( name );
    }
    return true;
  }
}

void PrimaryGeneratorAction::GeneratePrimaries( G4Event* event ) {
  if ( reader_ ) {
    if ( !generate_from_hepmc(event) ) {
      G4cout << "HepMC3 input exhausted after " << n_read_ << " input events;"
        " ending the run" << G4endl;
      G4RunManager::GetRunManager()->AbortRun( true );
    }
    return;
  }
  if ( isotropic_ ) {
    double cos_t = 2.*G4UniformRand() - 1.;
    double sin_t = std::sqrt( 1. - cos_t*cos_t );
    double phi = CLHEP::twopi * G4UniformRand();
    gun_->SetParticleMomentumDirection( G4ThreeVector(sin_t*std::cos(phi),
      sin_t*std::sin(phi), cos_t) );
  }
  gun_->GeneratePrimaryVertex( event );
}

void RunAction::EndOfRunAction( const G4Run* run ) {
  rec_->end_of_run( run );
}

void EventAction::BeginOfEventAction( const G4Event* event ) {
  rec_->begin_event( event );
}

void EventAction::EndOfEventAction( const G4Event* event ) {
  rec_->end_event( event );
}

void TrackingAction::PreUserTrackingAction( const G4Track* track ) {
  rec_->pre_track( track );
}

void TrackingAction::PostUserTrackingAction( const G4Track* track ) {
  rec_->post_track( track );
}

void SteppingAction::UserSteppingAction( const G4Step* step ) {
  rec_->step( step );
}
