/// @file
/// @brief Run-time options for the MARLEY de-excitation plug-in for Geant4
//
// This file is part of the MARLEY Geant4 de-excitation interface (g4marley).
// It is distributed under the same terms as MARLEY (GPL v3).

#pragma once

#include <cstdint>
#include <string>

/// How the initial spin-parity of a Geant4 fragment is chosen before it is
/// handed to MARLEY. Geant4 models (Bertini, PreCompound, BIC, ...) do not
/// track the nuclear spin of the residual nucleus, so some assumption is
/// needed. See README.md ("Physics choices") for a discussion.
enum class G4MarleySpinModel {
  /// Do not attach a spin-parity. MARLEY's NucleusDecayer then assigns one
  /// itself (in the GENIE-de-excite branch this is always J = 0 or 1/2 with
  /// positive parity, because the level-density weights are never filled).
  MarleyDefault,
  /// Sample (J, Pi) with probability proportional to MARLEY's own level
  /// density rho(Ex, J, Pi) for the fragment (the usual statistical-model
  /// assumption when only the excitation energy is known).
  LevelDensity,
  /// Always use the fixed values G4MarleyDeexConfig::fixed_twoJ / fixed_parity
  Fixed,
  /// If the fragment is exactly target + projectile (a compound nucleus that
  /// has emitted nothing yet, e.g. a few-MeV neutron absorbed by Bertini) and
  /// the projectile is n, p, d, t, 3He or alpha, sample (J, Pi) from the
  /// entrance channel: sigma(J,Pi) ~ (2J+1) sum_{l,j} T_lj(E_cm), with T_lj
  /// from MARLEY's optical model. Otherwise use LevelDensity.
  Auto
};

struct G4MarleyDeexConfig {

  /// Fragments with excitation energy below this value (MeV) are passed to
  /// the native Geant4 evaporation instead. MARLEY itself stops a continuum
  /// cascade at 1 keV, so going lower than that is pointless.
  double min_Ex = 1e-3;

  /// Fragments with Ex/A above this value (MeV per nucleon) are passed to the
  /// native Geant4 evaporation instead (MARLEY has no multifragmentation).
  double max_Ex_per_A = 10.;

  /// Fragments lighter than this are passed to the native Geant4 chain.
  /// Note that G4ExcitationHandler applies Fermi break-up to A < 17, Z < 9
  /// before evaporation is ever called, so those rarely reach MARLEY anyway.
  int min_A = 5;

  G4MarleySpinModel spin_model = G4MarleySpinModel::Auto;
  int fixed_twoJ = 0;
  int fixed_parity = 1;

  /// Largest 2J considered when sampling from the level density
  int max_twoJ = 30;

  /// Number of times a failed MARLEY decay (e.g., zero total width for the
  /// sampled J^pi, or not enough energy after the mass-table correction) is
  /// retried with a fresh spin-parity before falling back to Geant4.
  int max_attempts = 5;

  /// Seed for MARLEY's own random number generator. Zero means "derive a
  /// seed from the Geant4 random engine when the engine is created", so the
  /// full simulation is reproducible from the Geant4 seed alone.
  std::uint64_t seed = 0;

  /// 0 = silent, 1 = one line per MARLEY call, 2 = also dump the cascade
  int verbose = 0;

  /// Keep per-call records (needed for the HepMC3 de-excitation annotations
  /// and for validation output). Cheap; leave on unless you need speed.
  bool record = true;

  /// If true, every Bertini (G4CascadeInterface) instance is switched to its
  /// G4PreCompoundModel-based de-excitation so that its residual nuclei reach
  /// G4ExcitationHandler (and hence MARLEY). Without this, Bertini uses its
  /// internal evaporator and MARLEY is never called for Bertini events.
  bool bertini_use_precompound = true;

  /// For a compound nucleus that has emitted nothing yet (entrance-channel
  /// spin used), resample the MARLEY decay until the outcome is not compound
  /// elastic (projectile + target ground state). The hadronic inelastic
  /// cross section that triggered the interaction does not contain compound
  /// elastic, and Bertini rejects such outcomes anyway by re-running the
  /// whole cascade (up to 20 times, then returning "no interaction" and
  /// leaving the track unchanged). Doing the conditioning here is the same
  /// physics but converges and is far cheaper.
  bool cn_reject_elastic = true;
  /// Same, for pure gamma decay of a neutron-induced compound nucleus
  /// (radiative capture). Geant4 handles capture with a separate nCapture
  /// process and cross section, which the neutron inelastic cross section
  /// excludes. Not applied to charged projectiles, whose reaction cross
  /// sections include (p,gamma) etc.
  bool cn_reject_neutron_capture = true;
  /// Maximum number of resamplings for the two options above. If no
  /// acceptable outcome is found (channel closed), the last one is returned.
  int cn_max_resample = 300;

  /// If false, the plug-in is installed but every fragment goes to the
  /// native Geant4 evaporation. Useful for A/B comparisons with identical
  /// model routing (e.g., Bertini still using PreCompound).
  bool use_marley = true;
};
