add_library(trident_library
    Logging.cpp
    TridentHits.cpp
    filter.cpp
    fit.cpp
    process_g4_npe.cpp
    root_editor.cpp)

# We need this directory, and users of our library will need it too
target_include_directories(trident_library PUBLIC ../include)

target_link_libraries(trident_library PUBLIC ROOT::RIO ROOT::Tree ROOT::Gpad ROOT::Physics)
target_link_libraries(trident_library PUBLIC spdlog::spdlog)

