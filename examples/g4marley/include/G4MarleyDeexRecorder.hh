/// @file
/// @brief Bookkeeping of MARLEY de-excitation calls made from inside Geant4
//
// This file is part of the MARLEY Geant4 de-excitation interface (g4marley).
// It is distributed under the same terms as MARLEY (GPL v3).

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "G4LorentzVector.hh"

/// One binary step of a MARLEY cascade, as read back from MARLEY's HepMC3
/// record (vertex status 22 = Hauser-Feshbach continuum decay, 23 = discrete
/// gamma transition between tabulated levels).
struct G4MarleyDeexStep {
  int type = 0;
  int parent_pdg = 0;
  double parent_Ex = 0.;   // MeV
  int parent_twoJ = -1;
  int parent_parity = 0;
  int emitted_pdg = 0;
  double emitted_KE = 0.;  // MeV, in the rest frame of the initial fragment
  int daughter_pdg = 0;
  double daughter_Ex = 0.; // MeV
  int daughter_twoJ = -1;
  int daughter_parity = 0;
  double time_ns = 0.;     // decay time relative to fragment creation
};

/// A final-state particle handed back to Geant4
struct G4MarleyDeexProduct {
  int pdg = 0;
  G4LorentzVector p4;      // MeV, in the frame Geant4 gave the fragment in
  double time_ns = 0.;     // relative to the fragment creation time
  int emitter_pdg = 0;     // nucleus that emitted it (0 for the residual)
};

/// Everything about one MARLEY de-excitation of one Geant4 fragment
struct G4MarleyDeexRecord {
  long call_id = 0;
  int Z = 0;
  int A = 0;
  double Ex = 0.;          // MeV, Geant4's excitation energy
  int twoJ = -1;           // -1 = left for MARLEY to choose
  int parity = 0;
  int spin_source = 0;     // 0 MARLEY default, 1 level density,
                           // 2 entrance channel, 3 fixed
  G4LorentzVector p4;      // fragment four-momentum (MeV)
  std::vector<G4MarleyDeexProduct> products; // residual nucleus is last
  std::vector<G4MarleyDeexStep> steps;
  int attempts = 0;
  double lambda = 1.;      // momentum scale applied for exact E conservation
  double ex_left = 0.;     // MeV of MARLEY excitation left in the residual
                           // (absorbed into kinetic energy)
  double cpu_ms = 0.;
  int cn_resamples = 0;    // compound-elastic/capture outcomes resampled
  int proj_pdg = 0;        // entrance channel (if known), else 0
  double proj_KE = 0.;     // MeV, lab
  int target_pdg = 0;

  /// Compact one-line JSON (no spaces) used for HepMC3 attributes and logs
  std::string to_json() const;
  /// Human-readable cascade, e.g. "41Ar*(7.580,1/2+) -n-> 40Ar(gs)"
  std::string chain_string() const;
};

/// Per-thread singleton collecting records and counters. Records accumulate
/// in a pending list that the user application drains (e.g., once per G4
/// step) with TakePending().
class G4MarleyDeexRecorder {

  public:

    static G4MarleyDeexRecorder& Instance();

    struct Counters {
      long marley_calls = 0;        // fragments de-excited by MARLEY
      long marley_failed = 0;       // MARLEY failed max_attempts times
      long failed_attempts = 0;     // individual failed attempts
      long native_low_Ex = 0;       // Ex < min_Ex
      long native_high_Ex = 0;      // Ex/A > max_Ex_per_A
      long native_light = 0;        // A < min_A
      long native_disabled = 0;     // use_marley = false
      long cn_resamples = 0;        // compound elastic/capture resampled
      long cn_resample_failed = 0;  // no non-elastic outcome found
      double cpu_ms = 0.;           // time spent inside MARLEY
    };

    void add( G4MarleyDeexRecord&& rec );
    std::vector< G4MarleyDeexRecord > take_pending();
    void clear_pending() { pending_.clear(); }
    size_t num_pending() const { return pending_.size(); }

    Counters& counters() { return counters_; }
    const Counters& counters() const { return counters_; }

    void print_summary( std::ostream& os ) const;

    long next_call_id() { return ++last_call_id_; }

  private:

    G4MarleyDeexRecorder() = default;

    std::vector< G4MarleyDeexRecord > pending_;
    Counters counters_;
    long last_call_id_ = 0;
};

/// Nuclide label such as "40Ar", "4He", "n", "p", "gamma" for a PDG code
std::string g4marley_pdg_label( int pdg );
