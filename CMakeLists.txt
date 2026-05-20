cmake_minimum_required(VERSION 3.16)
project(reve_event_display CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ROOT must be built with: -Dwebgui=ON -Droot7=ON -Dgdml=ON
# (cvmfs LCG releases since LCG_103 already have these on).
#
# Component notes:
#   Core, RIO, Tree, MathCore - basic ROOT I/O and math.
#   Geom, Gdml                - TGeo geometry + GDML import.
#   ROOTEve                   - the ROOT7 EVE-7 library: REveManager,
#                               REveViewer, REveScene, REveGeoShape, and
#                               (if ever needed) the REveCalo* / projection
#                               classes.  The orthographic / perspective
#                               multi-viewer setup used by sea_cucumber
#                               needs only REveViewer, which lives here -- no
#                               extra component is required for it.
#   ROOTWebDisplay            - the web/local client transport EVE-7 uses to
#                               render viewers.
find_package(ROOT 6.30 CONFIG REQUIRED
    COMPONENTS Core Geom Gdml Tree RIO ROOTEve ROOTWebDisplay MathCore)

add_executable(event_display      event_display.cxx)
add_executable(sea_cucumber sea_cucumber.cxx)

set(_targets event_display sea_cucumber)
foreach(_t ${_targets})
    target_link_libraries(${_t} PRIVATE
        ROOT::Core
        ROOT::Geom
        ROOT::Gdml
        ROOT::Tree
        ROOT::RIO
        ROOT::ROOTEve
        ROOT::ROOTWebDisplay
        ROOT::MathCore
    )
    target_compile_options(${_t} PRIVATE
        -Wall -Wextra -Wpedantic
        $<$<CONFIG:Release>:-O2>
    )
endforeach()
