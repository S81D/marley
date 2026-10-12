/// @file
/// @brief G4VEvaporation implementation that hands nuclei to MARLEY
//
// This file is part of the MARLEY Geant4 de-excitation interface (g4marley).
// It is distributed under the same terms as MARLEY (GPL v3).

#pragma once

#include <memory>

#include "G4VEvaporation.hh"

class G4Evaporation;
class G4MarleyDeexEngine;

/// Drop-in replacement for G4Evaporation inside a G4ExcitationHandler.
///
/// G4ExcitationHandler::BreakItUp() calls BreakFragment() once per hot
/// fragment after (optional) multifragmentation and Fermi break-up. The
/// contract is: emitted products are appended to the vector, and the input
/// fragment is modified in place to become the residual nucleus. MARLEY
/// de-excites the nucleus all the way to its ground state in one call, so
/// the handler sees a cold residual and stops. Fragments MARLEY should not
/// handle (see G4MarleyDeexEngine::check) go to an internal G4Evaporation,
/// which is exactly what Geant4 would have used.
///
/// Ownership: the internal G4Evaporation owns the photon-evaporation channel
/// (as in stock Geant4). This class only mirrors the pointer so that
/// G4ExcitationHandler::SetEvaporation() picks up the same object, and nulls
/// it again before the base-class destructor would delete it.
class G4MarleyEvaporation : public G4VEvaporation {

  public:

    explicit G4MarleyEvaporation( std::shared_ptr< G4MarleyDeexEngine > engine );
    ~G4MarleyEvaporation() override;

    void BreakFragment( G4FragmentVector* products, G4Fragment* nucleus )
      override;

    void InitialiseChannels() override;

    void SetPhotonEvaporation( G4VEvaporationChannel* ptr ) override;

  private:

    void sync_fallback();

    std::shared_ptr< G4MarleyDeexEngine > engine_;
    G4Evaporation* fallback_;
};
