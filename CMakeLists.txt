cmake_minimum_required(VERSION 3.16)
project(reve_event_display CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ROOT must be built with: -Dwebgui=ON -Droot7=ON -Dgdml=ON
# (cvmfs LCG releases since LCG_103 already have these on).
find_package(ROOT 6.30 CONFIG REQUIRED
    COMPONENTS Core Geom Gdml Tree RIO ROOTEve ROOTWebDisplay)

add_executable(event_display event_display.cxx)

target_link_libraries(event_display PRIVATE
    ROOT::Core
    ROOT::Geom
    ROOT::Gdml
    ROOT::Tree
    ROOT::RIO
    ROOT::ROOTEve
    ROOT::ROOTWebDisplay
)

# Useful warnings
target_compile_options(event_display PRIVATE
    -Wall -Wextra -Wpedantic
    $<$<CONFIG:Release>:-O2>
)
