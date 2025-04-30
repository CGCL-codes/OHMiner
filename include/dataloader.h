#pragma once
#include<vector>
#include<map>
#include "hypergraph.h"
#include "pattern.h"

class Dataloader{
public:
    void load_DataGraph(hypergraph* &g, const char *hyperedge_data, const char *node_lable_data);
    void load_PatternGraph(Pattern* &p, const char *pattern_data, const char *pattern_node_lable);
    
};

