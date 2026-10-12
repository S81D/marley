/// @file
/// @brief Command-line options for the g4marley test application

#pragma once

#include <string>

#include "G4MarleyDeexConfig.hh"

struct AppOptions {
  // Run control
  int num_events = 10;
  long g4_seed = 12345;
  std::string macro;               // optional Geant4 macro (run before BeamOn)
  int g4_verbose = 0;

  // Physics
  std::string physics_list = "QGSP_BERT_HP"; // or BERT_ONLY (see main)
  std::string deex_mode = "marley";          // marley | native | off
  G4MarleyDeexConfig deex;

  // Geometry: cube of liquid argon centred on the origin
  double world_half_m = 5.;

  // Primaries from a HepMC3 file instead of the gun (e.g. GENIE final
  // states converted with scripts/gst_to_hepmc3.py). Status-1 particles of
  // each input event are injected from the gun position at t = 0.
  std::string hepmc_in_file;
  std::string hepmc_in_species = "nucleons"; // nucleons | no-leptons | all

  // Particle gun (defaults; /gun/ commands in a macro override them)
  std::string gun_particle = "neutron";
  double gun_energy_MeV = 14.;
  bool gun_isotropic = false;      // otherwise along +z

  // Output
  std::string hepmc_file = "g4marley.hepmc3"; // "none" disables
  std::string deex_log_file;       // JSON-lines log of every MARLEY call
  std::string had_csv_file;        // one row per secondary of every
                                   // non-elastic hadronic interaction
  double record_min_ekin_MeV = 0.; // e+/e-/gamma below this are not written
  bool store_trajectories = true;  // polyline points for charged tracks
};

/// Parses argv into opts. Returns false (after printing usage) on error.
bool parse_app_options( int argc, char* argv[], AppOptions& opts );
