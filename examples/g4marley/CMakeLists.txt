# g4marley: MARLEY nuclear de-excitation inside Geant4
#
# Builds
#   libG4MarleyDeex  - the plug-in (G4MarleyEvaporation + G4MarleyDeexPhysics)
#   g4marley         - test application with HepMC3 history output
#   g4marley_deex_test - feeds fixed fragments to G4ExcitationHandler (A/B)
#
# Needs Geant4 (10.6 tested), MARLEY (g4-deexcitation branch) and HepMC3.
# Point CMake at them with Geant4_DIR, MARLEY_DIR (or CMAKE_PREFIX_PATH) and
# HepMC3_DIR as usual. With no MARLEY_DIR, a working marley-config on PATH
# (from MARLEY's setup_marley.sh) is used first.

cmake_minimum_required(VERSION 3.16)
project(g4marley CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Geant4 REQUIRED)
find_package(HepMC3 REQUIRED)

# MARLEY can be found in two ways:
#  1. marley-config, which setup_marley.sh puts on PATH (and sets $MARLEY so
#     it reports the right directories). This is MARLEY's own mechanism for a
#     source-tree build and carries its include, library, GSL and ROOT flags.
#     Used unless MARLEY_DIR is given, provided its --libdir really holds
#     libMARLEY.
#  2. MARLEY's CMake package config (MARLEY_DIR or CMAKE_PREFIX_PATH). Two
#     things make find_package(MARLEY) fail hard, so they are checked first:
#      - The exported target names GSL:: and (if MARLEY was built with ROOT)
#        ROOT:: targets without looking for those packages itself, so we find
#        them here and skip the package config if CMake cannot.
#      - A build tree can hold a MARLEYConfig.cmake with no MARLEYTargets.cmake
#        beside it; that file is incomplete and is skipped.
find_program(MARLEY_CONFIG marley-config)
set(_marley_use_config_script FALSE)
if(MARLEY_CONFIG AND NOT MARLEY_DIR)
  execute_process(COMMAND ${MARLEY_CONFIG} --libdir
    OUTPUT_VARIABLE _marley_libdir OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(EXISTS "${_marley_libdir}/libMARLEY.so"
      OR EXISTS "${_marley_libdir}/libMARLEY.dylib")
    set(_marley_use_config_script TRUE)
  else()
    message(STATUS "Ignoring ${MARLEY_CONFIG}: no libMARLEY in its --libdir "
      "'${_marley_libdir}' (is $MARLEY set to the MARLEY source directory?)")
  endif()
endif()

if(NOT _marley_use_config_script)
  set(_marley_cfg_dir "")
  string(REPLACE ":" ";" _env_prefixes "$ENV{CMAKE_PREFIX_PATH}")
  foreach(_p IN LISTS MARLEY_DIR CMAKE_PREFIX_PATH _env_prefixes)
    if(_p STREQUAL "" OR _marley_cfg_dir)
      continue()
    endif()
    foreach(_sub "." "lib/cmake/MARLEY" "lib64/cmake/MARLEY")
      if(NOT _marley_cfg_dir AND EXISTS "${_p}/${_sub}/MARLEYConfig.cmake")
        if(EXISTS "${_p}/${_sub}/MARLEYTargets.cmake")
          get_filename_component(_marley_cfg_dir "${_p}/${_sub}" ABSOLUTE)
        else()
          message(STATUS "Ignoring incomplete MARLEY package config in "
            "${_p}/${_sub} (no MARLEYTargets.cmake)")
        endif()
      endif()
    endforeach()
  endforeach()

  if(_marley_cfg_dir)
    # Find the packages behind the GSL:: and ROOT:: targets MARLEY links to
    file(READ "${_marley_cfg_dir}/MARLEYTargets.cmake" _marley_targets)
    set(_marley_missing "")
    if(_marley_targets MATCHES "GSL::")
      find_package(GSL QUIET)
      if(NOT GSL_FOUND)
        list(APPEND _marley_missing GSL)
      endif()
    endif()
    string(REGEX MATCHALL "ROOT::[A-Za-z0-9_]+" _marley_root_targets
      "${_marley_targets}")
    if(_marley_root_targets)
      list(REMOVE_DUPLICATES _marley_root_targets)
      # ROOT::Core;ROOT::RIO -> Core;RIO
      string(REPLACE "ROOT::" "" _marley_root_comps "${_marley_root_targets}")
      find_package(ROOT QUIET COMPONENTS ${_marley_root_comps})
      if(NOT ROOT_FOUND)
        list(APPEND _marley_missing "ROOT (${_marley_root_comps})")
      endif()
    endif()
    if(_marley_missing)
      message(STATUS "Not using MARLEY package config in ${_marley_cfg_dir}: "
        "CMake cannot find ${_marley_missing}")
    else()
      find_package(MARLEY CONFIG QUIET PATHS "${_marley_cfg_dir}"
        NO_DEFAULT_PATH)
    endif()
  endif()
endif()

if(MARLEY_FOUND)
  message(STATUS "Found MARLEY via CMake config: ${MARLEY_DIR}")
elseif(NOT MARLEY_CONFIG)
  message(FATAL_ERROR "MARLEY not found. Source MARLEY's setup_marley.sh "
    "(puts marley-config on PATH) or pass -DMARLEY_DIR=<dir containing "
    "MARLEYConfig.cmake and MARLEYTargets.cmake>.")
else()
  execute_process(COMMAND ${MARLEY_CONFIG} --cflags
    OUTPUT_VARIABLE MARLEY_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${MARLEY_CONFIG} --libs
    OUTPUT_VARIABLE MARLEY_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
  separate_arguments(MARLEY_CFLAGS)
  separate_arguments(MARLEY_LIBS)
  # CMake adds no build RPATH for plain -L/-l flags (it does for the imported
  # MARLEY::MARLEY target), so add libMARLEY's directory by hand: the
  # binaries then run without LD_LIBRARY_PATH pointing at it.
  execute_process(COMMAND ${MARLEY_CONFIG} --libdir
    OUTPUT_VARIABLE _marley_libdir OUTPUT_STRIP_TRAILING_WHITESPACE)
  list(APPEND CMAKE_BUILD_RPATH "${_marley_libdir}")
  message(STATUS "Found MARLEY via marley-config: ${MARLEY_CONFIG}")
endif()

# ---- plug-in library ----
add_library(G4MarleyDeex SHARED
  src/G4MarleyDeexEngine.cc
  src/G4MarleyDeexRecorder.cc
  src/G4MarleyEvaporation.cc
  src/G4MarleyDeexPhysics.cc)
target_include_directories(G4MarleyDeex PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/include ${Geant4_INCLUDE_DIRS})
target_compile_definitions(G4MarleyDeex PUBLIC ${Geant4_DEFINITIONS})
target_link_libraries(G4MarleyDeex PUBLIC ${Geant4_LIBRARIES} HepMC3::HepMC3)
if(MARLEY_FOUND)
  target_link_libraries(G4MarleyDeex PUBLIC MARLEY::MARLEY)
else()
  target_compile_options(G4MarleyDeex PUBLIC ${MARLEY_CFLAGS})
  target_link_libraries(G4MarleyDeex PUBLIC ${MARLEY_LIBS})
endif()

# ---- test application ----
add_executable(g4marley
  app/g4marley.cc app/AppOptions.cc app/AppActions.cc app/HistoryRecorder.cc)
target_include_directories(g4marley PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/app)
target_link_libraries(g4marley PRIVATE G4MarleyDeex)

# ---- standalone de-excitation A/B test ----
add_executable(g4marley_deex_test tools/g4marley_deex_test.cc)
target_link_libraries(g4marley_deex_test PRIVATE G4MarleyDeex)

install(TARGETS G4MarleyDeex g4marley g4marley_deex_test
  LIBRARY DESTINATION lib RUNTIME DESTINATION bin)
install(DIRECTORY include/ DESTINATION include/g4marley)
