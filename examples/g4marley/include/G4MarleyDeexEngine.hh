/// @file
/// @brief Thin wrapper around marley::NucleusDecayer for Geant4 fragments
//
// This file is part of the MARLEY Geant4 de-excitation interface (g4marley).
// It is distributed under the same terms as MARLEY (GPL v3).

#pragma once

#include <memory>
#include <set>
#include <tuple>
#include <vector>

#include "G4LorentzVector.hh"

#include "G4MarleyDeexConfig.hh"
#include "G4MarleyDeexRecorder.hh"

class G4Fragment;

namespace marley {
  class Generator;
}

/// Owns a de-excitation-only marley::Generator and converts between Geant4
/// fragments and MARLEY's HepMC3 event record.
///
/// Frames and energy bookkeeping:
///  - MARLEY decays the nucleus at rest using *its* mass table, starting from
///    the excitation energy Geant4 assigned (G4Fragment::GetExcitationEnergy).
///  - The rest-frame three-momenta are then rescaled by a common factor
///    lambda so that the products, evaluated with *Geant4* masses, add up
///    exactly to the fragment's invariant mass. (The two codes' nuclear masses
///    differ by the electron binding energy, ~15 keV for argon, and by
///    mass-evaluation differences.) Momentum is conserved by construction
///    because the rest-frame momenta sum to zero.
///  - Finally everything is boosted back into the frame of the input fragment.
class G4MarleyDeexEngine {

  public:

    explicit G4MarleyDeexEngine( const G4MarleyDeexConfig& config );
    ~G4MarleyDeexEngine();

    const G4MarleyDeexConfig& config() const { return config_; }

    /// Why a fragment was not given to MARLEY (None means it should be)
    enum class Skip { None, LowEx, HighEx, Light, Disabled };
    Skip check( const G4Fragment& frag ) const;

    /// De-excite a fragment. On success fills the record (products in the
    /// fragment's frame, residual last) and returns true. On failure returns
    /// false and the caller should fall back to the native Geant4 model.
    bool decay( const G4Fragment& frag, G4MarleyDeexRecord& rec );

  private:

    bool try_decay_once( const G4Fragment& frag, G4MarleyDeexRecord& rec );

    /// Returns the spin source used (see G4MarleyDeexRecord::spin_source)
    int choose_spin_parity( const G4Fragment& frag, int& twoJ, int& parity );
    void level_density_spin_parity( int pdg, double Ex, int& twoJ,
      int& parity );
    /// Entrance-channel spin-parity for an un-decayed compound nucleus.
    /// Returns false if the current Geant4 interaction is not of that kind.
    bool entrance_channel_spin_parity( const G4Fragment& frag, int& twoJ,
      int& parity );

    /// True if a compound-nucleus decay produced an outcome that the
    /// triggering cross section excludes (see cn_reject_* options)
    bool unwanted_cn_outcome( const G4MarleyDeexRecord& rec ) const;

    /// One attempt-loop (up to max_attempts) of try_decay_once
    bool decay_with_attempts( const G4Fragment& frag, G4MarleyDeexRecord& rec,
      std::string& last_error );

    G4MarleyDeexConfig config_;
    std::unique_ptr< marley::Generator > gen_;

    /// Entrance channel found by the last entrance_channel_spin_parity call
    struct Entrance { int proj_pdg = 0; double proj_KE = 0.; int target_pdg = 0; };
    Entrance entrance_;

    /// Entrance channels (target, projectile, KE in keV) for which no
    /// acceptable compound-nucleus outcome was found. Bertini re-calls us up
    /// to 20 times with the same input, so a closed channel is only searched
    /// thoroughly once.
    std::set< std::tuple<int, int, long> > closed_channels_;
};
