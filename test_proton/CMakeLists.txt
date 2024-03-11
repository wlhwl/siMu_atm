cmake_minimum_required(VERSION 3.21)
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

find_package (corsika CONFIG REQUIRED)

# Find ROOT
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree)

# Include ROOT's headers
include(${ROOT_USE_FILE})

add_executable (c8_air_shower c8_air_shower.cpp initial_energy_generator.hpp)
target_link_libraries (c8_air_shower CORSIKA8::CORSIKA8)
target_link_libraries(c8_air_shower ${ROOT_LIBRARIES})

install (
  TARGETS c8_air_shower DESTINATION bin
)
