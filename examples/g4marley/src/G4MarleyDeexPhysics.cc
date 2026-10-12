/// @file
/// @brief Physics constructor that installs MARLEY de-excitation in Geant4
//
// This file is part of the MARLEY Geant4 de-excitation interface (g4marley).
// It is distributed under the same terms as MARLEY (GPL v3).

#include <set>

#include "G4CascadeInterface.hh"
#include "G4ExcitationHandler.hh"
#include "G4HadronicInteractionRegistry.hh"
#include "G4VPreCompoundModel.hh"
#include "G4ios.hh"

#include "G4MarleyDeexEngine.hh"
#include "G4MarleyDeexPhysics.hh"
#include "G4MarleyEvaporation.hh"

G4MarleyDeexPhysics::G4MarleyDeexPhysics( const G4MarleyDeexConfig& config,
  G4int verbose ) : G4VPhysicsConstructor( "MarleyDeexcitation" ),
  config_( config ), verbose_( verbose )
{
}

G4MarleyDeexPhysics::~G4MarleyDeexPhysics() = default;

void G4MarleyDeexPhysics::ConstructProcess() {

  auto* registry = G4HadronicInteractionRegistry::Instance();

  // ---- 1. Route Bertini's residual nuclei through G4ExcitationHandler ----
  int n_bertini = 0;
  if ( config_.bertini_use_precompound ) {
    for ( auto* model : registry->FindAllModels("BertiniCascade") ) {
      auto* casc = dynamic_cast< G4CascadeInterface* >( model );
      if ( !casc ) continue;
      // Creates a G4PreCompoundDeexcitation, which looks up (or creates and
      // registers) the shared "PRECO" model.
      casc->usePreCompoundDeexcitation();
      ++n_bertini;
    }
  }

  // ---- 2. Install MARLEY in every PRECO model's excitation handler ----
  // One engine (one marley::Generator) shared by all handlers on a thread.
  engine_ = std::make_shared< G4MarleyDeexEngine >( config_ );

  std::set< G4ExcitationHandler* > done;
  for ( auto* model : registry->FindAllModels("PRECO") ) {
    auto* preco = dynamic_cast< G4VPreCompoundModel* >( model );
    if ( !preco ) continue;
    G4ExcitationHandler* handler = preco->GetExcitationHandler();
    if ( !handler || done.count(handler) ) continue;
    // isLocal = true: the handler deletes our evaporation object, which in
    // turn deletes its internal G4Evaporation and photon evaporation.
    handler->SetEvaporation( new G4MarleyEvaporation(engine_), true );
    done.insert( handler );
  }

  if ( verbose_ > 0 ) {
    G4cout << "[g4marley] Bertini instances switched to PreCompound"
      << " de-excitation: " << n_bertini << '\n'
      << "[g4marley] G4ExcitationHandlers using MARLEY: " << done.size()
      << ( config_.use_marley ? "" : " (MARLEY disabled: native fallback only)" )
      << G4endl;
  }
  if ( done.empty() ) {
    G4cout << "[g4marley] WARNING: no PRECO model found; MARLEY will never"
      << " be called with this physics list" << G4endl;
  }
}
