cmake_minimum_required(VERSION 3.21)
project(analysis_atm_muon)

set(CMAKE_CXX_STANDARD 17)

file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/*.hpp)

#find ROOT
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Graf)

if(ROOT_FOUND)
    include(${ROOT_USE_FILE})
else()
    # Specify the ROOT headers directory
    include_directories(/home/ineffablord/documents/ROOT/root-6.28/include)
endif()

add_executable(analysis ${PROJECT_BINARY_DIR}/analysis.cc ${sources} ${headers})

if(ROOT_FOUND)
    target_link_libraries(analysis ${ROOT_LIBRARIES})
else()
    # Specify the ROOT libraries directory
    link_directories(/home/ineffablord/documents/ROOT/root-6.28/lib)
endif()