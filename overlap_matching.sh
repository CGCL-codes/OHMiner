#!/bin/bash

cd build
make -j
cd ..

dataset=("TC")
datasetDir="./hypergraph/hypergraph_unlable/"
patternDir="./pattern/"

EdgeNum=(2 3 4 5)
for data in "${dataset[@]}"
do
    for ((k=0;k<${#EdgeNum[*]};k++))
    do
        for ((i=1; i<=3; i++))
        do
            ./build/bin/overlap_match_parallel ${datasetDir}${data}"/hypergraph.txt" ${datasetDir}${data}"/node.txt"\
            ${patternDir}${data}"/patternEdge"${i}_${EdgeNum[k]} ${patternDir}${data}"/patternNode"${i}_${EdgeNum[k]}
        done
    done
done