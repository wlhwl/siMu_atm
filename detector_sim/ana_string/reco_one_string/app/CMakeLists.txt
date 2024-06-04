add_executable(analysis analysis.cpp)
#target_compile_features(apps PRIVATE cxx_std_17)

target_link_libraries(analysis PRIVATE trident_library)

