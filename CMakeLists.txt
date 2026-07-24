# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(Geant4-MicroLight-BNZS)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build  with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${PROJECT_SOURCE_DIR}/include/core
  ${PROJECT_SOURCE_DIR}/include/modes
  ${PROJECT_SOURCE_DIR}/include/stageA
  ${PROJECT_SOURCE_DIR}/include/stageB
  ${PROJECT_SOURCE_DIR}/include/stageD
  )


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB_RECURSE sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB_RECURSE headers ${PROJECT_SOURCE_DIR}/include/*.hh)

# Keep the v3 RVE reader explicit: it is a required geometry component rather
# than an optional source discovered only through the recursive glob.
list(APPEND sources ${PROJECT_SOURCE_DIR}/src/core/RvePlacementReader.cc)
list(APPEND headers ${PROJECT_SOURCE_DIR}/include/core/RvePlacementReader.hh)
list(REMOVE_DUPLICATES sources)
list(REMOVE_DUPLICATES headers)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
set(BNZS_EXECUTABLE Geant4-MicroLight-BNZS)
add_executable(${BNZS_EXECUTABLE} Geant4MicroLightBNZS.cc ${sources} ${headers})
target_link_libraries(${BNZS_EXECUTABLE} ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build the executable. This is so that we can run it directly because it
# relies on these scripts being in the current working directory.
#
set(BNZS_SCRIPTS
  init_vis.mac
  vis.mac
  StageA.mac
  StageD_OpticalHomogenization.mac
  )

foreach(_script ${BNZS_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
#  standalone
#
#add_custom_target(${BNZS_EXECUTABLE} DEPENDS exampleB1)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS ${BNZS_EXECUTABLE} DESTINATION bin)
