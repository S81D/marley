/// @file
/// @brief Bookkeeping of MARLEY de-excitation calls made from inside Geant4
//
// This file is part of the MARLEY Geant4 de-excitation interface (g4marley).
// It is distributed under the same terms as MARLEY (GPL v3).

#include <cstdio>
#include <sstream>

#include "G4MarleyDeexRecorder.hh"

namespace {

  const char* const ELEMENT_SYMBOLS[] = { "n", "H", "He", "Li", "Be", "B",
    "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
    "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La",
    "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
    "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb",
    "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U" };
  constexpr int NUM_SYMBOLS = sizeof( ELEMENT_SYMBOLS ) / sizeof( const char* );

  // Formats a number compactly for JSON (no spaces, enough digits for keV
  // resolution on MeV-scale quantities)
  std::string num( double x ) {
    char buf[32];
    std::snprintf( buf, sizeof(buf), "%.8g", x );
    return buf;
  }

  std::string spin_parity( int twoJ, int parity ) {
    if ( twoJ < 0 ) return "";
    std::string s = ( twoJ % 2 == 0 ) ? std::to_string( twoJ / 2 )
      : std::to_string( twoJ ) + "/2";
    if ( parity > 0 ) s += "+";
    else if ( parity < 0 ) s += "-";
    return s;
  }

  std::string state_label( int pdg, double Ex, int twoJ, int parity ) {
    std::string s = g4marley_pdg_label( pdg );
    char buf[32];
    if ( Ex > 1e-6 ) {
      std::snprintf( buf, sizeof(buf), "%.3f", Ex );
      s += "*(" + std::string( buf );
      std::string jp = spin_parity( twoJ, parity );
      if ( !jp.empty() ) s += "," + jp;
      s += ")";
    }
    else {
      std::string jp = spin_parity( twoJ, parity );
      s += "(gs" + ( jp.empty() ? std::string() : "," + jp ) + ")";
    }
    return s;
  }

}

std::string g4marley_pdg_label( int pdg ) {
  switch ( pdg ) {
    case 22: return "gamma";
    case 11: return "e-";
    case -11: return "e+";
    case 2112: return "n";
    case 2212: return "p";
    case 1000010020: return "d";
    case 1000010030: return "t";
    case 1000020030: return "3He";
    case 1000020040: return "alpha";
    default: break;
  }
  if ( pdg > 1000000000 ) {
    int Z = ( pdg / 10000 ) % 1000;
    int A = ( pdg / 10 ) % 1000;
    std::string sym = ( Z < NUM_SYMBOLS ) ? ELEMENT_SYMBOLS[Z]
      : "Z" + std::to_string( Z );
    return std::to_string( A ) + sym;
  }
  return std::to_string( pdg );
}

std::string G4MarleyDeexRecord::chain_string() const {
  if ( steps.empty() ) return "";
  std::ostringstream os;
  os << state_label( steps.front().parent_pdg, steps.front().parent_Ex,
    steps.front().parent_twoJ, steps.front().parent_parity );
  for ( const auto& s : steps ) {
    char buf[32];
    std::snprintf( buf, sizeof(buf), "%.3f", s.emitted_KE );
    os << " -" << g4marley_pdg_label( s.emitted_pdg ) << "(" << buf
      << ")-> " << state_label( s.daughter_pdg, s.daughter_Ex,
      s.daughter_twoJ, s.daughter_parity );
  }
  return os.str();
}

std::string G4MarleyDeexRecord::to_json() const {
  std::ostringstream os;
  os << "{\"id\":" << call_id << ",\"Z\":" << Z << ",\"A\":" << A
    << ",\"Ex\":" << num( Ex ) << ",\"twoJ\":" << twoJ << ",\"par\":"
    << parity << ",\"jsrc\":" << spin_source << ",\"p4\":[" << num( p4.px() ) << ',' << num( p4.py() ) << ','
    << num( p4.pz() ) << ',' << num( p4.e() ) << "],\"lam\":"
    << num( lambda ) << ",\"exl\":" << num( ex_left ) << ",\"att\":"
    << attempts << ",\"ms\":" << num( cpu_ms ) << ",\"ncn\":" << cn_resamples
    << ",\"proj\":" << proj_pdg << ",\"projKE\":" << num( proj_KE )
    << ",\"tgt\":" << target_pdg << ",\"prod\":[";
  for ( size_t i = 0; i < products.size(); ++i ) {
    const auto& p = products[i];
    if ( i ) os << ',';
    os << '[' << p.pdg << ',' << num( p.p4.px() ) << ',' << num( p.p4.py() )
      << ',' << num( p.p4.pz() ) << ',' << num( p.p4.e() ) << ','
      << num( p.time_ns ) << ',' << p.emitter_pdg << ']';
  }
  os << "],\"steps\":[";
  for ( size_t i = 0; i < steps.size(); ++i ) {
    const auto& s = steps[i];
    if ( i ) os << ',';
    os << '[' << s.type << ',' << s.parent_pdg << ',' << num( s.parent_Ex )
      << ',' << s.parent_twoJ << ',' << s.parent_parity << ','
      << s.emitted_pdg << ',' << num( s.emitted_KE ) << ',' << s.daughter_pdg
      << ',' << num( s.daughter_Ex ) << ',' << s.daughter_twoJ << ','
      << s.daughter_parity << ',' << num( s.time_ns ) << ']';
  }
  os << "]}";
  return os.str();
}

G4MarleyDeexRecorder& G4MarleyDeexRecorder::Instance() {
  // One recorder per thread, matching one engine per thread
  static thread_local G4MarleyDeexRecorder instance;
  return instance;
}

void G4MarleyDeexRecorder::add( G4MarleyDeexRecord&& rec ) {
  pending_.push_back( std::move(rec) );
}

std::vector< G4MarleyDeexRecord > G4MarleyDeexRecorder::take_pending() {
  std::vector< G4MarleyDeexRecord > out;
  out.swap( pending_ );
  return out;
}

void G4MarleyDeexRecorder::print_summary( std::ostream& os ) const {
  const auto& c = counters_;
  os << "==== g4marley de-excitation summary ====\n"
     << "  fragments de-excited by MARLEY : " << c.marley_calls << '\n'
     << "  MARLEY failures (-> native G4) : " << c.marley_failed
     << "  (failed attempts: " << c.failed_attempts << ")\n"
     << "  native G4, Ex < min_Ex          : " << c.native_low_Ex << '\n'
     << "  native G4, Ex/A > max           : " << c.native_high_Ex << '\n'
     << "  native G4, A < min_A            : " << c.native_light << '\n'
     << "  native G4, MARLEY disabled      : " << c.native_disabled << '\n'
     << "  compound-nucleus resamples      : " << c.cn_resamples
     << "  (no acceptable outcome: " << c.cn_resample_failed << ")\n"
     << "  CPU time inside MARLEY          : " << c.cpu_ms / 1000. << " s";
  if ( c.marley_calls > 0 ) {
    os << " (" << c.cpu_ms / c.marley_calls << " ms/call)";
  }
  os << "\n========================================\n";
}
