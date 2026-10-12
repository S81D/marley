/// @file
/// @brief Geometry, primary generator and user actions for g4marley

#pragma once

#include <map>
#include <memory>
#include <string>

#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "AppOptions.hh"

class G4ParticleGun;
class G4ParticleDefinition;
class HistoryRecorder;
namespace HepMC3 { class ReaderAscii; }

/// Information about the input (generator) event behind the current Geant4
/// event, filled by PrimaryGeneratorAction and written by HistoryRecorder
struct GenEventInfo {
  long input_index = -1;                       // position in the input file
  int input_event_number = -1;
  std::map< std::string, std::string > attrs;  // input event attributes
};
GenEventInfo& CurrentGenEventInfo();

/// A cube of liquid argon (NIST G4_lAr, natural isotopic composition)
class DetectorConstruction : public G4VUserDetectorConstruction {
  public:
    explicit DetectorConstruction( double half_size_m ) : half_m_( half_size_m ) {}
    G4VPhysicalVolume* Construct() override;
  private:
    double half_m_;
};

/// Particle gun at the origin, along +z or isotropic (/gun/ commands work),
/// or primaries read from a HepMC3 file (--hepmc-in).
class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  public:
    explicit PrimaryGeneratorAction( const AppOptions& opts );
    ~PrimaryGeneratorAction() override;
    void GeneratePrimaries( G4Event* event ) override;
  private:
    bool generate_from_hepmc( G4Event* event );
    bool keep_species( int pdg ) const;
    const G4ParticleDefinition* definition( int pdg ) const;

    G4ParticleGun* gun_;
    bool isotropic_;
    std::unique_ptr< HepMC3::ReaderAscii > reader_;
    std::string species_;
    long n_read_ = 0;
};

class RunAction : public G4UserRunAction {
  public:
    explicit RunAction( std::shared_ptr< HistoryRecorder > rec ) : rec_( rec ) {}
    void EndOfRunAction( const G4Run* run ) override;
  private:
    std::shared_ptr< HistoryRecorder > rec_;
};

class EventAction : public G4UserEventAction {
  public:
    explicit EventAction( std::shared_ptr< HistoryRecorder > rec ) : rec_( rec ) {}
    void BeginOfEventAction( const G4Event* event ) override;
    void EndOfEventAction( const G4Event* event ) override;
  private:
    std::shared_ptr< HistoryRecorder > rec_;
};

class TrackingAction : public G4UserTrackingAction {
  public:
    explicit TrackingAction( std::shared_ptr< HistoryRecorder > rec )
      : rec_( rec ) {}
    void PreUserTrackingAction( const G4Track* track ) override;
    void PostUserTrackingAction( const G4Track* track ) override;
  private:
    std::shared_ptr< HistoryRecorder > rec_;
};

class SteppingAction : public G4UserSteppingAction {
  public:
    explicit SteppingAction( std::shared_ptr< HistoryRecorder > rec )
      : rec_( rec ) {}
    void UserSteppingAction( const G4Step* step ) override;
  private:
    std::shared_ptr< HistoryRecorder > rec_;
};
