#include "common.h"
#include "hypergraph.h"
#include "pattern.h"
#include "dataloader.h"
#include "algorithm"
#include "iostream"
#include <immintrin.h>
#include <x86intrin.h>
double get_wall_time() {
    struct timeval time;
    if(gettimeofday(&time,NULL)) {
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * 0.000001;
}

int binarySerach(const std::vector<int>&array, int key) {
    int left = 0;
    int right = array.size() - 1;
    while (left <= right) {
        int mid = (left + right) / 2;
        if (array[mid] == key) {
            return mid;
        }
        else if (array[mid] < key) {
            left = mid + 1;
        }
        else {
            right = mid - 1;
        }
    }
    return -1;
}


bool has_data(const std::vector<int>&set ,int depth, int v){
    for(int i = 0; i < depth; ++i){
        if(set[i] == v)
            return true;
    }
    return false;
}

std::vector<int> sub_merge(std::vector<int>&a, std::vector<int>&b){
    std::vector<int>re;
    int i = 0, j = 0;
    int size1 = a.size(), size2 = b.size();
    while(i < size1 && j < size2){
        int v1 = a[i], v2 = b[j];
        if(v1 > v2){
            ++j;
        }
        else if(v1 < v2){
            re.push_back(a[i]);
            ++i;
        }
        else{
            ++j;
            ++i;
        }
    }
    while(i < size1)
        re.push_back(a[i++]);
    return re;
}

std::vector<int> intersect(const std::vector<int>&set1, const std::vector<int>&set2){
    std::vector<int>ans;
    auto size1 = set1.size(), size2 = set2.size();
    ans.reserve(std::min(size1, size2));
    if(size1 != 0 && size2 != 0 && size1 / size2 > 10){
        for(int i = 0; i < size2; ++i){
            int index = binarySerach(set1, set2[i]);
            if (index == -1)
                continue;
            else
                ans.emplace_back(set2[i]);
        }
        return ans;
    }
    else if(size1 != 0 && size2 != 0 && size2 / size1 > 10){
        for(int i = 0; i < size1; ++i){
            int index = binarySerach(set2, set1[i]);
            if (index == -1)
                continue;
            else
                ans.emplace_back(set1[i]);
        }
        return ans;
    }
    auto i = 0, j = 0;
    while(i < size1 && j < size2){
        if(set1[i] < set2[j])
            ++i;
        else if(set1[i] > set2[j])
            ++j;
        else{
            ans.push_back(set1[i]);
            ++i;
            ++j;
        }
    }
    return ans;
}

std::vector<int> intersect(const std::vector<int>&set1, const std::vector<int>&set2, int& l, int& r){
    std::vector<int>ans;
    auto i = 0, j = l;
    auto size1 = set1.size();
    auto size2 = r;
    while(i < size1 && j < size2){
        if(set1[i] < set2[j])
            ++i;
        else if(set1[i] > set2[j])
            ++j;
        else{
            ans.push_back(set1[i]);
            ++i;
            ++j;
        }
    }
    return ans;
}



bool intersect_overlap(const std::vector<int>&set1, const std::vector<int>&set2){
    int i = 0, j = 0, size1 = set1.size(), size2 = set2.size();
    while(i < size1 && j < size2){
        if(set1[i] < set2[j])
            ++i;
        else if(set1[i] > set2[j])
            ++j;
        else{
            return true;
        }
    }
    return false;
}



std::vector<int> unordered_sub_hash(const std::vector<int>& set1, const std::vector<int>& set2){
    auto size1 = set1.size();
    std::vector<int> re;
    for (int j = 0; j < size1; ++j)
        if (!std::binary_search(set2.begin(), set2.end() , set1[j]))
            re.emplace_back(set1[j]);
    return re;
}

void unordered_sub_hash(const std::vector<int>& set1, const std::vector<int>& set2, std::vector<int>&re){
    auto size1 = set1.size();
    
    for (int j = 0; j < size1; ++j)
        if (!std::binary_search(set2.begin(), set2.end() , set1[j]))
            re.emplace_back(set1[j]);
    return;
}

std::vector<int> unordered_sub_hash(const std::vector<int>& set1, const std::vector<int>& set2, int& l, int& r){
    auto size1 = set1.size();
    std::vector<int> re;
    for (int j = 0; j < size1; ++j)
        if (!std::binary_search(set2.begin() + l, set2.begin() + r, set1[j]))
            re.emplace_back(set1[j]);
    return re;
}


//多个集合求交集
void intersection(const std::vector<std::set<int>>& candidate_union, std::vector<int>& candidate) {
    if (candidate_union.empty()) {
        return;
    }
    candidate = std::vector<int>(candidate_union.begin()->begin(), candidate_union.begin()->end());

    // 依次计算集合的交集
    for (auto it = std::next(candidate_union.begin()); it != candidate_union.end(); ++it) {
        std::set<int> temp;
        std::set_intersection(candidate.begin(), candidate.end(),
                              it->begin(), it->end(),
                              std::inserter(temp, temp.begin()));
        candidate = std::vector<int>(temp.begin(), temp.end());
    }
}

void intersect(const std::vector<int>&set1, const std::vector<int>&set2, std::vector<int>& set3){
    set3.clear();
    auto size1 = set1.size(), size2 = set2.size();
    set3.reserve(std::min(size1, size2));
    if(size1 != 0 && size2 != 0 && size1 / size2 > 10){
        for(int i = 0; i < size2; ++i){
            int index = binarySerach(set1, set2[i]);
            if (index == -1)
                continue;
            else
                set3.emplace_back(set2[i]);
        }
        return;
    }
    else if(size1 != 0 && size2 != 0 && size2 / size1 > 10){
        for(int i = 0; i < size1; ++i){
            int index = binarySerach(set2, set1[i]);
            if (index == -1)
                continue;
            else
                set3.emplace_back(set1[i]);
        }
        return;
    }
    auto i = 0, j = 0;
    while(i < size1 && j < size2){
        if(set1[i] < set2[j])
            ++i;
        else if(set1[i] > set2[j])
            ++j;
        else{
            set3.emplace_back(set1[i]);
            ++i;
            ++j;
        }
    }
    return;
}


