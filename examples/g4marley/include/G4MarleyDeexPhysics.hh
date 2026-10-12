/// @file
/// @brief Physics constructor that installs MARLEY de-excitation in Geant4
//
// This file is part of the MARLEY Geant4 de-excitation interface (g4marley).
// It is distributed under the same terms as MARLEY (GPL v3).

#pragma once

#include <memory>

#include "G4VPhysicsConstructor.hh"

#include "G4MarleyDeexConfig.hh"

class G4MarleyDeexEngine;

/// Register this *after* the reference physics list has been created, e.g.
///
///   G4VModularPhysicsList* list = factory.GetReferencePhysList("QGSP_BERT_HP");
///   list->RegisterPhysics( new G4MarleyDeexPhysics(config) );
///
/// Constructors run ConstructProcess() in registration order, so by the time
/// ours runs every hadronic model already exists but no G4ExcitationHandler
/// has been initialised yet (that happens when physics tables are built at
/// the first BeamOn). We then
///   1. optionally switch every Bertini instance to PreCompound
///      de-excitation (so its residual nuclei reach G4ExcitationHandler), and
///   2. give every shared "PRECO" model's G4ExcitationHandler a
///      G4MarleyEvaporation.
/// No Geant4 source modification is needed.
///
/// Models covered: anything that de-excites through the registry's "PRECO"
/// G4PreCompoundModel: Bertini (with step 1), BIC, INCL++, FTF/QGS string
/// models (via G4GeneratorPrecompoundInterface), mu- capture at rest
/// (G4MuMinusCapturePrecompound), and PreCompound used directly.
/// NOT covered: NeutronHP/ParticleHP (it uses evaluated final states and its
/// own photon handling), G4NeutronRadCapture, QMD, models that build a
/// private G4ExcitationHandler, and Bertini if step 1 is disabled.
class G4MarleyDeexPhysics : public G4VPhysicsConstructor {

  public:

    explicit G4MarleyDeexPhysics( const G4MarleyDeexConfig& config,
      G4int verbose = 1 );
    ~G4MarleyDeexPhysics() override;

    void ConstructParticle() override {}
    void ConstructProcess() override;

    const std::shared_ptr< G4MarleyDeexEngine >& engine() const
      { return engine_; }

  private:

    G4MarleyDeexConfig config_;
    G4int verbose_;
    std::shared_ptr< G4MarleyDeexEngine > engine_;
};
