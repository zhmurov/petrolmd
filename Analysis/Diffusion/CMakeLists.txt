add_executable(diffusion main.cpp)

target_compile_features(diffusion PRIVATE cxx_std_20)

target_link_libraries(diffusion PUBLIC libchemfiles.a mdio)