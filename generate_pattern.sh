#!/bin/bash

dataset=("CP" "HB" "SB" "WT" "TC")
datasetDir="/dataset/graph/hypergraph_unlable/"
patternDir="/dataset/pattern/pattern_unlable/"


restrict=("2 5 15" "3 10 20")

EdgeNum=(2)

cd build
make -j
cd ..

for data in "${dataset[@]}"
do
    mkdir -p ${patternDir}${data}
    for ((k=0;k<${#EdgeNum[*]};k++))
    do
        for ((i=1; i<=3; i++))
        do
            ./build/bin/random_sample ${datasetDir}${data}"/hypergraph.txt" ${datasetDir}${data}"/node.txt" ${patternDir}${data}"/patternEdge"${i}_${EdgeNum[k]} ${patternDir}${data}"/patternNode"${i}_${EdgeNum[k]} ${restrict[k]}
            sleep 2
        done
    done
done