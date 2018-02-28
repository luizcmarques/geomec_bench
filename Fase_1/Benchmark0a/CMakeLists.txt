file(GLOB headers *.h)
file(GLOB sources *.cpp)

set(Mymain main.cpp ${sources} ${headers})

add_executable(Benchmark0a ${Mymain})

target_link_libraries(Benchmark0a pz)