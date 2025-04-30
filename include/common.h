#pragma once
#include <sys/time.h>
#include <cstdlib>
#include <string>
#include "vector"
#include "set"

double get_wall_time();
int binarySerach(const std::vector<int>&array, int key);
bool has_data(const std::vector<int>&set ,int depth, int v);
std::vector<int> sub_merge(std::vector<int>&a, std::vector<int>&b);
std::vector<int> intersect(const std::vector<int>&set1, const std::vector<int>&set2);
void intersect(const std::vector<int>&set1, const std::vector<int>&set2, std::vector<int>& set3);
std::vector<int> intersect(const std::vector<int>&set1, const std::vector<int>&set2, int& l, int& r);
bool intersect_overlap(const std::vector<int>&set1, const std::vector<int>&set2);
std::vector<int> unordered_sub_hash(const std::vector<int>& set1, const std::vector<int>& set2);
void unordered_sub_hash(const std::vector<int>& set1, const std::vector<int>& set2, std::vector<int>&re);
std::vector<int> unordered_sub_hash(const std::vector<int>& set1, const std::vector<int>& set2, int& l, int& r);
void intersection(const std::vector<std::set<int>>& candidate_union, std::vector<int>& candidate);
std::vector<int> intersect_simd(const std::vector<int>& larray, const std::vector<int>& rarray);
std::vector<int> intersect_simd512(const std::vector<int>& larray, const std::vector<int>& rarray);
void intersect_simd512(const std::vector<int>& larray, const std::vector<int>& rarray, std::vector<int>&cn);