/// @file
/// @brief Command-line options for the g4marley test application

#include <iostream>
#include <stdexcept>
#include <string>

#include "AppOptions.hh"

namespace {

  void print_usage( const char* exe ) {
    std::cout << "Usage: " << exe << " [options] [macro.mac]\n"
      "  -n N                 number of events (default 10)\n"
      "  -o FILE              HepMC3 output file, or 'none' (default g4marley.hepmc3)\n"
      "  --deex-log FILE      JSON-lines log of every MARLEY call\n"
      "  --had-csv FILE       CSV, one row per secondary of every non-elastic\n"
      "                       hadronic interaction (compact; for big samples)\n"
      "  -p LIST              reference physics list (default QGSP_BERT_HP);\n"
      "                       BERT_ONLY = EM + Bertini neutron inelastic only\n"
      "                       (reproduces the original gpvm test setup)\n"
      "  --deex MODE          marley | native | off (default marley)\n"
      "                         native: plug-in installed, all fragments to G4\n"
      "                         off:    stock Geant4, plug-in not installed\n"
      "  --bertini-deex MODE  precompound | cascade (default precompound)\n"
      "  --spin MODEL         auto | ld | marley | fixed:2J:P (default auto)\n"
      "                         auto: entrance channel for an un-decayed compound\n"
      "                         nucleus from n,p,d,t,3He,alpha; else level density\n"
      "  --min-ex MEV         minimum Ex handed to MARLEY (default 0.001)\n"
      "  --no-cn-conditioning let MARLEY return compound-elastic/capture outcomes\n"
      "                       for an un-decayed compound nucleus (Bertini will\n"
      "                       then retry the whole cascade, as in the gpvm log)\n"
      "  --max-ex-per-a MEV   maximum Ex/A handed to MARLEY (default 10)\n"
      "  --seed N             Geant4 seed (MARLEY seed derived from it)\n"
      "  --marley-seed N      explicit MARLEY seed\n"
      "  --world-half-m X     half-size of the LAr cube in m (default 5)\n"
      "  --hepmc-in FILE      take primaries from a HepMC3 file (status-1 particles\n"
      "                       of each event, injected at the gun position, t = 0)\n"
      "  --hepmc-in-species S nucleons | no-leptons | all (default nucleons)\n"
      "  --gun-particle NAME  (default neutron)\n"
      "  --gun-energy MEV     (default 14)\n"
      "  --gun-iso            isotropic gun direction (default +z)\n"
      "  --record-min-ekin E  don't write e+/e-/gamma with KE < E MeV\n"
      "  --no-traj            don't store trajectory points\n"
      "  -v N                 g4marley verbosity (1 = one line per MARLEY call)\n"
      "  --g4verbose N        Geant4 tracking/run verbosity\n";
  }

}

bool parse_app_options( int argc, char* argv[], AppOptions& o ) {
  try {
    for ( int i = 1; i < argc; ++i ) {
      std::string a = argv[i];
      auto next = [&]() -> std::string {
        if ( i + 1 >= argc ) throw std::runtime_error( "missing value for "
          + a );
        return argv[++i];
      };
      if ( a == "-h" || a == "--help" ) { print_usage( argv[0] ); return false; }
      else if ( a == "-n" ) o.num_events = std::stoi( next() );
      else if ( a == "-o" ) o.hepmc_file = next();
      else if ( a == "--deex-log" ) o.deex_log_file = next();
      else if ( a == "--had-csv" ) o.had_csv_file = next();
      else if ( a == "-p" ) o.physics_list = next();
      else if ( a == "--deex" ) {
        o.deex_mode = next();
        if ( o.deex_mode != "marley" && o.deex_mode != "native"
          && o.deex_mode != "off" ) throw std::runtime_error( "bad --deex" );
      }
      else if ( a == "--bertini-deex" ) {
        std::string m = next();
        if ( m == "precompound" ) o.deex.bertini_use_precompound = true;
        else if ( m == "cascade" ) o.deex.bertini_use_precompound = false;
        else throw std::runtime_error( "bad --bertini-deex" );
      }
      else if ( a == "--spin" ) {
        std::string m = next();
        if ( m == "auto" ) o.deex.spin_model = G4MarleySpinModel::Auto;
        else if ( m == "ld" ) o.deex.spin_model = G4MarleySpinModel::LevelDensity;
        else if ( m == "marley" ) {
          o.deex.spin_model = G4MarleySpinModel::MarleyDefault;
        }
        else if ( m.rfind("fixed:", 0) == 0 ) {
          // fixed:2J:P, e.g. fixed:1:1 for J^pi = 1/2+
          size_t c2 = m.find( ':', 6 );
          if ( c2 == std::string::npos ) throw std::runtime_error( "bad --spin" );
          o.deex.spin_model = G4MarleySpinModel::Fixed;
          o.deex.fixed_twoJ = std::stoi( m.substr(6, c2 - 6) );
          o.deex.fixed_parity = std::stoi( m.substr(c2 + 1) ) >= 0 ? 1 : -1;
        }
        else throw std::runtime_error( "bad --spin" );
      }
      else if ( a == "--min-ex" ) o.deex.min_Ex = std::stod( next() );
      else if ( a == "--no-cn-conditioning" ) {
        o.deex.cn_reject_elastic = false;
        o.deex.cn_reject_neutron_capture = false;
      }
      else if ( a == "--max-ex-per-a" ) o.deex.max_Ex_per_A = std::stod( next() );
      else if ( a == "--seed" ) o.g4_seed = std::stol( next() );
      else if ( a == "--marley-seed" ) o.deex.seed = std::stoull( next() );
      else if ( a == "--world-half-m" ) o.world_half_m = std::stod( next() );
      else if ( a == "--hepmc-in" ) o.hepmc_in_file = next();
      else if ( a == "--hepmc-in-species" ) {
        o.hepmc_in_species = next();
        if ( o.hepmc_in_species != "nucleons" && o.hepmc_in_species != "no-leptons"
          && o.hepmc_in_species != "all" ) throw std::runtime_error(
          "bad --hepmc-in-species" );
      }
      else if ( a == "--gun-particle" ) o.gun_particle = next();
      else if ( a == "--gun-energy" ) o.gun_energy_MeV = std::stod( next() );
      else if ( a == "--gun-iso" ) o.gun_isotropic = true;
      else if ( a == "--record-min-ekin" ) {
        o.record_min_ekin_MeV = std::stod( next() );
      }
      else if ( a == "--no-traj" ) o.store_trajectories = false;
      else if ( a == "-v" ) o.deex.verbose = std::stoi( next() );
      else if ( a == "--g4verbose" ) o.g4_verbose = std::stoi( next() );
      else if ( !a.empty() && a[0] != '-' ) o.macro = a;
      else throw std::runtime_error( "unknown option " + a );
    }
  }
  catch ( const std::exception& e ) {
    std::cerr << "Error: " << e.what() << "\n";
    print_usage( argv[0] );
    return false;
  }
  o.deex.use_marley = ( o.deex_mode == "marley" );
  return true;
}
