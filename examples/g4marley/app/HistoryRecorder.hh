/// @file
/// @brief Writes the full Geant4 particle history of each event to HepMC3

#pragma once

#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "G4ThreeVector.hh"
#include "G4VUserTrackInformation.hh"

#include "HepMC3/GenEvent.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/WriterAscii.h"

#include "AppOptions.hh"
#include "G4MarleyDeexRecorder.hh"

class G4Event;
class G4Run;
class G4Step;
class G4Track;

/// Links a secondary G4Track (before it has a track ID) to the HepMC3
/// particle created for it in the parent's step
class HistoryTrackInfo : public G4VUserTrackInformation {
  public:
    HistoryTrackInfo( HepMC3::GenParticlePtr p, bool rec )
      : particle( std::move(p) ), recorded( rec ) {}
    HepMC3::GenParticlePtr particle;
    bool recorded;
};

/// Event record conventions (also written into the HepMC3 run info):
///  * One HepMC3 particle per *segment* of a Geant4 track. A track is split
///    into a new segment (same "tid" attribute) wherever it has an
///    interaction vertex but survives (elastic scatter, Compton, delta ray..)
///  * Vertex = point where a track interacted / produced secondaries.
///    Position (x, y, z, c*t) in mm. Attribute "proc" = Geant4 process name.
///    Vertex status: 11 hadronic, 12 decay, 13 EM, 14 other. Primaries
///    hang off the event root vertex (event position = gun position).
///  * Particle status: 4 primary (first segment), 2 ended in a vertex,
///    1 no end vertex (left the world, stopped, or below tracking cuts).
///  * Particle attributes: tid, ptid (parent tid), cproc (creator process),
///    eproc (process that ended a track without an end vertex,
///    "Transportation" = left the world),
///    traj (x,y,z,ct quadruplets in mm: intermediate points of charged tracks
///    and the end point of tracks without an end vertex), deex (MARLEY call
///    id if the particle is a MARLEY de-excitation product).
///  * Vertex attribute "marley": JSON array with the MARLEY cascade(s)
///    (see G4MarleyDeexRecord::to_json) for de-excitations at that vertex.
class HistoryRecorder {

  public:

    HistoryRecorder( const AppOptions& opts );
    ~HistoryRecorder();

    void begin_event( const G4Event* event );
    void end_event( const G4Event* event );
    void pre_track( const G4Track* track );
    void post_track( const G4Track* track );
    void step( const G4Step* step );
    void end_of_run( const G4Run* run );

    /// Extra run-level metadata (set before the first event)
    void set_run_attribute( const std::string& key, const std::string& value )
      { run_attrs_[key] = value; }

  private:

    struct OpenTrack {
      HepMC3::GenParticlePtr particle;
      bool recorded = false;
      std::vector< double > traj;
      G4ThreeVector last_pos;
      G4ThreeVector last_dir;
    };

    HepMC3::GenParticlePtr make_particle( const G4Track* trk, int status ) const;
    bool keep_secondary( const G4Track* sec ) const;
    void push_point( OpenTrack& o, const G4ThreeVector& pos, double t_ns );
    void finish_segment( OpenTrack& o, bool has_end_vertex );

    /// Matches MARLEY records made during this step to the step's
    /// secondaries. Returns the accepted records; matched secondaries get
    /// their call id in matched_ids (parallel to the secondaries vector).
    std::vector< const G4MarleyDeexRecord* > match_marley(
      const std::vector< G4MarleyDeexRecord >& recs,
      const std::vector< const G4Track* >& secs,
      std::vector< long >& matched_ids ) const;

    void write_deex_log( const G4Step* step,
      const std::vector< G4MarleyDeexRecord >& recs,
      const std::vector< const G4MarleyDeexRecord* >& accepted );

    AppOptions opts_;
    std::map< std::string, std::string > run_attrs_;
    std::shared_ptr< HepMC3::GenRunInfo > run_info_;
    std::unique_ptr< HepMC3::WriterAscii > writer_;
    std::unique_ptr< HepMC3::GenEvent > event_;
    std::unordered_map< int, OpenTrack > open_;
    std::ofstream deex_log_;
    std::ofstream had_csv_;
    long had_vertex_counter_ = 0;
    int current_event_ = -1;

    // Per-event and per-run diagnostics
    long evt_marley_accepted_ = 0;
    long evt_marley_discarded_ = 0;
    long run_marley_accepted_ = 0;
    long run_marley_discarded_ = 0;
    long run_no_interaction_steps_ = 0;  // hadronic step, nothing happened
    std::map< int, long > calls_per_step_; // histogram over hadronic steps
};
