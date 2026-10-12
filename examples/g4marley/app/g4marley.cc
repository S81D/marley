/// @file
/// @brief g4marley: Geant4 test application with MARLEY nuclear de-excitation
///
/// Shoots particles into a cube of liquid argon, de-excites residual nuclei
/// with MARLEY (via G4MarleyDeexPhysics), and writes the full particle
/// history of every event to a HepMC3 ASCII file. Run with -h for options.

#include <iostream>
#include <memory>
#include <string>

#include "G4CascadeInterface.hh"
#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4HadronicProcess.hh"
#include "G4Neutron.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4PhysListFactory.hh"
#include "G4PhysicsListHelper.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VPhysicsConstructor.hh"
#include "Randomize.hh"

#include "AppActions.hh"
#include "AppOptions.hh"
#include "G4MarleyDeexEngine.hh"
#include "G4MarleyDeexPhysics.hh"
#include "HistoryRecorder.hh"

namespace {

  /// Neutron inelastic scattering with Bertini from 0 eV, and nothing else
  /// hadronic (no elastic, no capture, no HP). This mirrors the hadronic
  /// summary printed in the original gpvm log (out.txt), and is only useful
  /// for reproducing that problem.
  class BertiniNeutronOnlyPhysics : public G4VPhysicsConstructor {
    public:
      BertiniNeutronOnlyPhysics() : G4VPhysicsConstructor( "BertNeutronOnly" ) {}
      void ConstructParticle() override {}
      void ConstructProcess() override {
        auto* proc = new G4NeutronInelasticProcess( "neutronInelastic" );
        proc->AddDataSet( new G4NeutronInelasticXS() );
        auto* bertini = new G4CascadeInterface();
        bertini->SetMinEnergy( 0. );
        bertini->SetMaxEnergy( 10.*GeV );
        proc->RegisterMe( bertini );
        G4PhysicsListHelper::GetPhysicsListHelper()->RegisterProcess( proc,
          G4Neutron::Neutron() );
      }
  };

  class BertiniOnlyList : public G4VModularPhysicsList {
    public:
      BertiniOnlyList() {
        RegisterPhysics( new G4EmStandardPhysics() );
        RegisterPhysics( new G4DecayPhysics() );
        RegisterPhysics( new BertiniNeutronOnlyPhysics() );
      }
  };

  const char* spin_model_name( G4MarleySpinModel m ) {
    switch ( m ) {
      case G4MarleySpinModel::MarleyDefault: return "marley-default";
      case G4MarleySpinModel::LevelDensity: return "level-density";
      case G4MarleySpinModel::Fixed: return "fixed";
      case G4MarleySpinModel::Auto: return "auto";
    }
    return "?";
  }

}

int main( int argc, char* argv[] ) {

  AppOptions opts;
  if ( !parse_app_options(argc, argv, opts) ) return 1;

  G4Random::setTheSeed( opts.g4_seed );

  auto* rm = new G4RunManager;
  rm->SetUserInitialization( new DetectorConstruction(opts.world_half_m) );

  G4VModularPhysicsList* list = nullptr;
  if ( opts.physics_list == "BERT_ONLY" ) list = new BertiniOnlyList();
  else {
    G4PhysListFactory factory;
    list = factory.GetReferencePhysList( opts.physics_list );
    if ( !list ) {
      std::cerr << "Unknown physics list " << opts.physics_list << '\n';
      return 2;
    }
  }

  // Install MARLEY de-excitation (must be registered after the hadronic
  // physics so that all models exist when its ConstructProcess() runs)
  G4MarleyDeexPhysics* deex_phys = nullptr;
  if ( opts.deex_mode != "off" ) {
    deex_phys = new G4MarleyDeexPhysics( opts.deex );
    list->RegisterPhysics( deex_phys );
  }
  rm->SetUserInitialization( list );

  auto recorder = std::make_shared< HistoryRecorder >( opts );
  rm->SetUserAction( new PrimaryGeneratorAction(opts) );
  rm->SetUserAction( new RunAction(recorder) );
  rm->SetUserAction( new EventAction(recorder) );
  rm->SetUserAction( new TrackingAction(recorder) );
  rm->SetUserAction( new SteppingAction(recorder) );

  rm->Initialize();

  auto* ui = G4UImanager::GetUIpointer();
  if ( opts.g4_verbose > 0 ) {
    ui->ApplyCommand( "/run/verbose " + std::to_string(opts.g4_verbose) );
    ui->ApplyCommand( "/tracking/verbose "
      + std::to_string(opts.g4_verbose - 1) );
  }
  if ( !opts.macro.empty() ) ui->ApplyCommand( "/control/execute "
    + opts.macro );

  recorder->set_run_attribute( "physics_list", opts.physics_list );
  recorder->set_run_attribute( "deex_mode", opts.deex_mode );
  recorder->set_run_attribute( "bertini_deex",
    opts.deex.bertini_use_precompound && opts.deex_mode != "off"
    ? "precompound" : "cascade" );
  recorder->set_run_attribute( "spin_model",
    spin_model_name(opts.deex.spin_model) );
  recorder->set_run_attribute( "g4_seed", std::to_string(opts.g4_seed) );
  if ( deex_phys && deex_phys->engine() ) {
    recorder->set_run_attribute( "marley_seed",
      std::to_string(deex_phys->engine()->config().seed) );
  }
  if ( opts.hepmc_in_file.empty() ) {
    recorder->set_run_attribute( "gun", opts.gun_particle + " "
      + std::to_string(opts.gun_energy_MeV) + " MeV"
      + ( opts.gun_isotropic ? " isotropic" : " +z" ) );
  }
  else {
    recorder->set_run_attribute( "primaries", "HepMC3 input " + opts.hepmc_in_file
      + ", species " + opts.hepmc_in_species + ", injected at the origin" );
  }
  recorder->set_run_attribute( "world", "G4_lAr cube, half-size "
    + std::to_string(opts.world_half_m) + " m" );

  std::cout << "Simulating " << opts.num_events << " events . . .\n";
  rm->BeamOn( opts.num_events );

  delete rm;
  return 0;
}
