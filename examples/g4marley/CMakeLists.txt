# g4marley: MARLEY nuclear de-excitation inside Geant4
#
# Builds
#   libG4MarleyDeex  - the plug-in (G4MarleyEvaporation + G4MarleyDeexPhysics)
#   g4marley         - test application with HepMC3 history output
#   g4marley_deex_test - feeds fixed fragments to G4ExcitationHandler (A/B)
#
# Needs Geant4 (10.6 tested), MARLEY (g4-deexcitation branch) and HepMC3.
# Point CMake at them with Geant4_DIR, MARLEY_DIR (or CMAKE_PREFIX_PATH) and
# HepMC3_DIR as usual.

cmake_minimum_required(VERSION 3.16)
project(g4marley CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Geant4 REQUIRED)
find_package(HepMC3 REQUIRED)

# MARLEY: prefer its CMake package config; fall back to marley-config.
# MARLEY's exported target refers to GSL::gsl but does not look for GSL
# itself, so look for GSL first, and only use the package config if GSL was
# found (otherwise configuring fails on the missing GSL::gsl target).
# marley-config carries its own GSL flags.
find_package(GSL QUIET)
if(GSL_FOUND)
  find_package(MARLEY QUIET)
else()
  message(STATUS "GSL not found by CMake: using marley-config for MARLEY")
endif()
if(MARLEY_FOUND)
  message(STATUS "Found MARLEY via CMake config: ${MARLEY_DIR}")
else()
  find_program(MARLEY_CONFIG marley-config REQUIRED)
  execute_process(COMMAND ${MARLEY_CONFIG} --cflags
    OUTPUT_VARIABLE MARLEY_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${MARLEY_CONFIG} --libs
    OUTPUT_VARIABLE MARLEY_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
  separate_arguments(MARLEY_CFLAGS)
  separate_arguments(MARLEY_LIBS)
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
