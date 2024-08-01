cmake_minimum_required(VERSION 3.21)
project(analysis_atm_muon)

set(CMAKE_CXX_STANDARD 17)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hpp)

#find ROOT
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Graf)

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

#git clone https://github.com/gabime/spdlog.git
#cd spdlog && mkdir build && cd build
#cmake -DCMAKE_INSTALL_PREFIX="$HOME/.local" ..
#make -j && make install
set(spdlog_DIR $ENV{HOME}/.local/lib/cmake/spdlog) # check if it is installed at $HOME/.local/lib or lib64...
find_package(spdlog REQUIRED)

if(ROOT_FOUND)
    include(${ROOT_USE_FILE})
else()
    # Specify the ROOT headers directory
    include_directories(/home/ineffablord/documents/ROOT/root-6.28/include)
endif()

#add_executable(analysis ${PROJECT_SOURCE_DIR}/app/analysis.cc ${sources} ${headers})

#if(ROOT_FOUND)
#    target_link_libraries(analysis ${ROOT_LIBRARIES})
#else()
#    # Specify the ROOT libraries directory
#    link_directories(/home/ineffablord/documents/ROOT/root-6.28/lib)
#endif()

add_subdirectory(src)
add_subdirectory(app)
