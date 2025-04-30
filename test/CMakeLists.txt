include_directories(../include)

add_executable(overlap_match overlap_match.cpp)
TARGET_LINK_LIBRARIES(overlap_match HyperMining "/opt/intel/oneapi/vtune/latest/lib64/libittnotify.a" -ldl)

add_executable(overlap_match_parallel overlap_match_parallel.cpp)
TARGET_LINK_LIBRARIES(overlap_match_parallel HyperMining "/opt/intel/oneapi/vtune/latest/lib64/libittnotify.a" -ldl)

add_executable(random_sample random_sample.cpp)
TARGET_LINK_LIBRARIES(random_sample HyperMining)