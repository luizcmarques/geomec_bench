file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Benchmark0a main.cpp ${headers})

add_executable(Benchmark0a ${Benchmark0a})

target_link_libraries(Benchmark0a pz)