/// @file
/// @brief Feeds identical excited nuclei to G4ExcitationHandler::BreakItUp()
///        with either MARLEY or native Geant4 evaporation installed.
///
/// Usage: g4marley_deex_test Z A Ex_MeV N MODE OUT_PREFIX [options]
///   MODE = marley | native
///   options: --spin ld|marley|fixed:2J:P   --seed N   --verbose N
///
/// Writes OUT_PREFIX_products.csv (one row per emitted particle) and
/// OUT_PREFIX_decays.csv (one row per decay with energy/momentum balance).
/// In marley mode also OUT_PREFIX_steps.csv (the MARLEY cascade steps).
/// No pre-equilibrium stage is involved: this is the equilibrium
/// de-excitation alone, which is exactly what the plug-in replaces.

#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>

#include "G4BaryonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4ExcitationHandler.hh"
#include "G4Fragment.hh"
#include "G4GenericIon.hh"
#include "G4IonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4NucleiProperties.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4ReactionProductVector.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include "G4MarleyDeexEngine.hh"
#include "G4MarleyDeexRecorder.hh"
#include "G4MarleyEvaporation.hh"

int main( int argc, char* argv[] ) {

  if ( argc < 7 ) {
    std::cout << "Usage: " << argv[0] << " Z A Ex_MeV N marley|native"
      " OUT_PREFIX [--spin ld|marley|fixed:2J:P] [--seed N] [--verbose N]\n";
    return 1;
  }
  const int Z = std::stoi( argv[1] );
  const int A = std::stoi( argv[2] );
  const double Ex = std::stod( argv[3] );
  const long N = std::stol( argv[4] );
  const std::string mode = argv[5];
  const std::string prefix = argv[6];

  G4MarleyDeexConfig cfg;
  long seed = 4242;
  for ( int i = 7; i < argc; ++i ) {
    std::string a = argv[i];
    if ( a == "--seed" && i + 1 < argc ) seed = std::stol( argv[++i] );
    else if ( a == "--verbose" && i + 1 < argc ) cfg.verbose = std::stoi( argv[++i] );
    else if ( a == "--spin" && i + 1 < argc ) {
      std::string m = argv[++i];
      if ( m == "auto" ) cfg.spin_model = G4MarleySpinModel::Auto;
      else if ( m == "ld" ) cfg.spin_model = G4MarleySpinModel::LevelDensity;
      else if ( m == "marley" ) cfg.spin_model = G4MarleySpinModel::MarleyDefault;
      else if ( m.rfind("fixed:", 0) == 0 ) {
        size_t c2 = m.find( ':', 6 );
        cfg.spin_model = G4MarleySpinModel::Fixed;
        cfg.fixed_twoJ = std::stoi( m.substr(6, c2 - 6) );
        cfg.fixed_parity = std::stoi( m.substr(c2 + 1) ) >= 0 ? 1 : -1;
      }
    }
  }
  G4Random::setTheSeed( seed );

  // Minimal particle setup so that G4IonTable can create ions
  G4BosonConstructor().ConstructParticle();
  G4LeptonConstructor().ConstructParticle();
  G4MesonConstructor().ConstructParticle();
  G4BaryonConstructor().ConstructParticle();
  G4IonConstructor().ConstructParticle();
  G4ShortLivedConstructor().ConstructParticle();
  G4GenericIon* gion = G4GenericIon::GenericIon();
  gion->SetProcessManager( new G4ProcessManager(gion) );
  G4ParticleTable::GetParticleTable()->SetReadiness();

  auto* handler = new G4ExcitationHandler();
  std::shared_ptr< G4MarleyDeexEngine > engine;
  if ( mode == "marley" ) {
    engine = std::make_shared< G4MarleyDeexEngine >( cfg );
    handler->SetEvaporation( new G4MarleyEvaporation(engine), true );
  }
  else if ( mode != "native" ) {
    std::cerr << "MODE must be marley or native\n";
    return 1;
  }

  const double M = G4NucleiProperties::GetNuclearMass( A, Z ) + Ex;
  const G4LorentzVector p4_initial( 0., 0., 0., M );

  std::ofstream prod_csv( prefix + "_products.csv" );
  std::ofstream dec_csv( prefix + "_decays.csv" );
  std::ofstream step_csv;
  prod_csv << "decay,pdg,Z,A,KE\n";
  dec_csv << "decay,n_products,dE,dp,residual_Z,residual_A,ms\n";
  if ( engine ) {
    step_csv.open( prefix + "_steps.csv" );
    step_csv << "decay,type,parent_pdg,parent_Ex,parent_twoJ,parent_parity,"
      "emitted_pdg,emitted_KE,daughter_pdg,daughter_Ex,daughter_twoJ,"
      "daughter_parity,time_ns\n";
  }

  auto& recorder = G4MarleyDeexRecorder::Instance();
  std::map< std::string, long > mult;

  for ( long d = 0; d < N; ++d ) {
    G4Fragment frag( A, Z, p4_initial );
    auto t0 = std::chrono::steady_clock::now();
    G4ReactionProductVector* products = handler->BreakItUp( frag );
    double ms = std::chrono::duration< double, std::milli >(
      std::chrono::steady_clock::now() - t0 ).count();

    G4LorentzVector sum;
    const size_t n_products = products->size();
    int rZ = 0, rA = 0, best_A = -1;
    for ( auto* rp : *products ) {
      const auto* def = rp->GetDefinition();
      G4LorentzVector lv( rp->GetMomentum(), rp->GetTotalEnergy() );
      sum += lv;
      int pa = def->GetAtomicMass();
      int pz = def->GetAtomicNumber();
      prod_csv << d << ',' << def->GetPDGEncoding() << ',' << pz << ',' << pa
        << ',' << rp->GetKineticEnergy()/MeV << '\n';
      if ( pa > best_A ) { best_A = pa; rZ = pz; rA = pa; }
      ++mult[ def->GetParticleName() ];
      delete rp;
    }
    delete products;
    dec_csv << d << ',' << n_products << ','
      << ( sum.e() - M )/MeV << ',' << sum.vect().mag()/MeV << ',' << rZ << ','
      << rA << ',' << ms << '\n';

    if ( engine ) {
      for ( const auto& rec : recorder.take_pending() ) {
        for ( const auto& s : rec.steps ) {
          step_csv << d << ',' << s.type << ',' << s.parent_pdg << ','
            << s.parent_Ex << ',' << s.parent_twoJ << ',' << s.parent_parity
            << ',' << s.emitted_pdg << ',' << s.emitted_KE << ','
            << s.daughter_pdg << ',' << s.daughter_Ex << ','
            << s.daughter_twoJ << ',' << s.daughter_parity << ','
            << s.time_ns << '\n';
        }
      }
    }
  }

  std::cout << "Z=" << Z << " A=" << A << " Ex=" << Ex << " MeV, " << N
    << " decays, mode " << mode << "\nAverage multiplicities:\n";
  for ( const auto& kv : mult ) {
    std::cout << "  " << kv.first << ": " << double( kv.second ) / N << '\n';
  }
  if ( engine ) recorder.print_summary( std::cout );
  delete handler;
  return 0;
}
