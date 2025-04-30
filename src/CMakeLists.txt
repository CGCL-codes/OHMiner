include_directories(../include)
INCLUDE_DIRECTORIES(${OpenMP_INCLUDE_DIRS})
set(HyperGraphSrc
        dataloader.cpp
        hypergraph.cpp
        common.cpp
        pattern.cpp)

add_library(HyperMining
            SHARED ${HyperGraphSrc})