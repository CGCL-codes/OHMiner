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

std::vector<int> intersect_simd(const std::vector<int>& larray, const std::vector<int>& rarray) {
    std::vector<int> cn;
    if (larray.empty() || rarray.empty())
        return cn;

    size_t l_count = larray.size();
    size_t r_count = rarray.size();

    // 如果larray的大小大于rarray的大小，则交换l_count和r_count的值
    if (l_count > r_count) {
        return intersect_simd(rarray, larray); // 交换形参并重新调用函数
    }
    cn.reserve(l_count);

    __m256i per_u_order = _mm256_set_epi32(1, 1, 1, 1, 0, 0, 0, 0);
    __m256i per_v_order = _mm256_set_epi32(3, 2, 1, 0, 3, 2, 1, 0);

    size_t li = 0;
    size_t ri = 0;


    auto size_ratio = (r_count) / (l_count);
    if (size_ratio > 2) {
        while (li < l_count && ri + 7 < r_count) {
            __m256i u_elements = _mm256_set1_epi32(larray[li]);
            __m256i v_elements = _mm256_loadu_si256((__m256i *) (rarray.data() + ri));

            while (true) {
                __m256i mask = _mm256_cmpeq_epi32(u_elements, v_elements);
                auto real_mask = _mm256_movemask_epi8(mask);
                if (real_mask != 0) {
                    cn.push_back(larray[li]);
                }
                if (larray[li] > rarray[ri + 7]) {
                    ri += 8;
                    if (ri + 7 >= r_count) {
                        break;
                    }
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray.data() + ri));
                } else {
                    li++;
                    if (li >= l_count) {
                        break;
                    }
                    u_elements = _mm256_set1_epi32(larray[li]);
                }
            }
        }
    } else {
        while (li + 1 < l_count && ri + 3 < r_count) {
            __m256i u_elements = _mm256_loadu_si256((__m256i *) (larray.data() + li));
            __m256i u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
            __m256i v_elements = _mm256_loadu_si256((__m256i *) (rarray.data() + ri));
            __m256i v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);

            while (true) {
                __m256i mask = _mm256_cmpeq_epi32(u_elements_per, v_elements_per);
                auto real_mask = _mm256_movemask_epi8(mask);
                if (real_mask << 16 != 0) {
                    cn.push_back(larray[li]);
                }
                if (real_mask >> 16 != 0) {
                    cn.push_back(larray[li + 1]);
                }

                if (larray[li + 1] == rarray[ri + 3]) {
                    li += 2;
                    ri += 4;
                    if (li + 1 >= l_count || ri + 3 >= r_count) {
                        break;
                    }
                    u_elements = _mm256_loadu_si256((__m256i *) (larray.data() + li));
                    u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray.data() + ri));
                    v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);
                } else if (larray[li + 1] > rarray[ri + 3]) {
                    ri += 4;
                    if (ri + 3 >= r_count) {
                        break;
                    }
                    v_elements = _mm256_loadu_si256((__m256i *) (rarray.data() + ri));
                    v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);
                } else {
                    li += 2;
                    if (li + 1 >= l_count) {
                        break;
                    }
                    u_elements = _mm256_loadu_si256((__m256i *) (larray.data() + li));
                    u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
                }
            }
        }
    }


    if (li < l_count && ri < r_count) {
        while (true) {
            while (larray[li] < rarray[ri]) {
                ++li;
                if (li >= l_count) {
                    return cn;
                }
            }
            while (larray[li] > rarray[ri]) {
                ++ri;
                if (ri >= r_count) {
                    return cn;
                }
            }
            if (larray[li] == rarray[ri]) {
                // 写回交集元素
                cn.push_back(larray[li]);

                ++li;
                ++ri;
                if (li >= l_count || ri >= r_count) {
                    return cn;
                }
            }
        }
    }

    return cn;
}

std::vector<int> intersect_simd512(const std::vector<int>& larray, const std::vector<int>& rarray) {
    int l_count = larray.size(), r_count = rarray.size();
    std::vector<int>cn;
    if (l_count == 0 || r_count == 0)
        return cn;

    int lc = l_count;
    int rc = r_count;

    if (lc > rc) {
        return intersect_simd512(rarray, larray);
    }
    cn.reserve(std::min(l_count, r_count));
    int li = 0;
    int ri = 0;

    // __m512i st = _mm512_set_epi32(3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 4, 3, 2, 1);
    __m512i st = _mm512_set_epi32(3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0);


    // __m512i p4f128 = _mm512_set_epi64(1, 0, 1, 0, 1, 0, 1, 0);
    __m512i p4f128 = _mm512_set_epi32(3,2,1, 0, 3,2,1, 0, 3,2,1, 0,3,2, 1, 0);
                    // auto x = _mm512_permutexvar_epi64(p4f128, st);// _mm512_permutex_epi64(v_elements, 0b01000100);;
                    // unsigned long long int a[8];
                    // _mm512_store_epi64(a, x);
                    // for (int i = 0; i< 8; ++i) {
                    //   printf("%llx\n", a[i]);
                    // }
                    //   return ;

    int* cur_back_ptr = cn.data();

    auto size1 = (rc) / (lc);
    if (size1 > 2) {
        if (li < lc && ri + 15 < rc) {
            __m512i u_elements = _mm512_set1_epi32(larray[li]);
            __m512i v_elements = _mm512_loadu_si512((__m512i *) (rarray.data() + ri));

            while (true) {
                __mmask16 mask = _mm512_cmpeq_epi32_mask(u_elements, v_elements);
                if (mask != 0x0000) {
                    // write back
                    _mm512_mask_compressstoreu_epi32(cur_back_ptr, mask, u_elements);
                    cur_back_ptr += _popcnt32(mask);
                }

                if (larray[li] > rarray[ri + 15]) {
                    ri += 16;
                    if (ri + 15 >= rc) {
                        break;
                    }
                    v_elements = _mm512_loadu_si512((__m512i *) (rarray.data() + ri));
                } else {
                    li += 1;
                    if (li >= lc) {
                        break;
                    }
                    u_elements = _mm512_set1_epi32(larray[li]);
                }
            }
        }
    } else {
        if (li + 3 < lc && ri + 3 < rc) {
            __m512i u_elements = _mm512_loadu_si512((__m512i *) (larray.data() + li));
            __m512i u_elements_per = _mm512_permutexvar_epi32(st, u_elements);
            __m512i v_elements = _mm512_loadu_si512((__m512i *) (rarray.data() + ri));
            // __m512i v_elements_per = _mm512_permutex_epi64(v_elements, 0b01000100);
            __m512i v_elements_per = _mm512_permutexvar_epi32(p4f128, v_elements);

            while (true) {
                __mmask16 mask = _mm512_cmpeq_epi32_mask(u_elements_per, v_elements_per);
                if (mask != 0x0000) {
                    // write back
                    _mm512_mask_compressstoreu_epi32(cur_back_ptr, mask, u_elements_per);
                    cur_back_ptr += _popcnt32(mask);
                }

                if (larray[li + 3] > rarray[ri + 3]) {
                    ri += 4;
                    if (ri + 3 >= rc) {
                        break;
                    }
                    v_elements = _mm512_loadu_si512((__m512i *) (rarray.data() + ri));
                    v_elements_per = _mm512_permutexvar_epi32(p4f128, v_elements);// _mm512_permutex_epi64(v_elements, 0b01000100);;
                } else if (larray[li + 3] < rarray[ri + 3]) {
                    li += 4;
                    if (li + 3 >= lc) {
                        break;
                    }
                    u_elements = _mm512_loadu_si512((__m512i *) (larray.data() + li));
                    u_elements_per = _mm512_permutexvar_epi32(st, u_elements);
                } else {
                    li += 4;
                    ri += 4;
                    if (li + 3 >= lc || ri + 3 >= rc) {
                        break;
                    }
                    u_elements = _mm512_loadu_si512((__m512i *) (larray.data() + li));
                    u_elements_per = _mm512_permutexvar_epi32(st, u_elements);
                    v_elements = _mm512_loadu_si512((__m512i *) (rarray.data() + ri));
                    v_elements_per = _mm512_permutexvar_epi32(p4f128, v_elements);//_mm512_permutex_epi64(v_elements, 0b01000100);
                }
            }
        }
    }

    cn.assign(cn.data(), cur_back_ptr);

    if (li < lc && ri < rc) {
        while (true) {
            while (larray[li] < rarray[ri]) {
                li += 1;
                if (li >= lc) {
                    return cn;
                }
            }
            while (larray[li] > rarray[ri]) {
                ri += 1;
                if (ri >= rc) {
                    return cn;
                }
            }
            if (larray[li] == rarray[ri]) {
                // write back
                //cn[cn_count++] = larray[li];
                cn.emplace_back(larray[li]);
                li += 1;
                ri += 1;
                if (li >= lc || ri >= rc) {
                    return cn;
                }
            }
        }
    }
    return cn;
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

void intersect_simd512(const std::vector<int>& larray, const std::vector<int>& rarray, std::vector<int>&cn) {
    int l_count = larray.size(), r_count = rarray.size();
    if (l_count == 0 || r_count == 0)
        return;

    int lc = l_count;
    int rc = r_count;

    if (lc > rc) {
        intersect_simd512(rarray, larray, cn);
        return;
    }
    cn.reserve(std::min(l_count, r_count));
    int li = 0;
    int ri = 0;

    // __m512i st = _mm512_set_epi32(3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 4, 3, 2, 1);
    __m512i st = _mm512_set_epi32(3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0);


    // __m512i p4f128 = _mm512_set_epi64(1, 0, 1, 0, 1, 0, 1, 0);
    __m512i p4f128 = _mm512_set_epi32(3,2,1, 0, 3,2,1, 0, 3,2,1, 0,3,2, 1, 0);
                    // auto x = _mm512_permutexvar_epi64(p4f128, st);// _mm512_permutex_epi64(v_elements, 0b01000100);;
                    // unsigned long long int a[8];
                    // _mm512_store_epi64(a, x);
                    // for (int i = 0; i< 8; ++i) {
                    //   printf("%llx\n", a[i]);
                    // }
                    //   return ;

    int* cur_back_ptr = cn.data();

    auto size1 = (rc) / (lc);
    if (size1 > 2) {
        if (li < lc && ri + 15 < rc) {
            __m512i u_elements = _mm512_set1_epi32(larray[li]);
            __m512i v_elements = _mm512_loadu_si512((__m512i *) (rarray.data() + ri));

            while (true) {
                __mmask16 mask = _mm512_cmpeq_epi32_mask(u_elements, v_elements);
                if (mask != 0x0000) {
                    // write back
                    _mm512_mask_compressstoreu_epi32(cur_back_ptr, mask, u_elements);
                    cur_back_ptr += _popcnt32(mask);
                }

                if (larray[li] > rarray[ri + 15]) {
                    ri += 16;
                    if (ri + 15 >= rc) {
                        break;
                    }
                    v_elements = _mm512_loadu_si512((__m512i *) (rarray.data() + ri));
                } else {
                    li += 1;
                    if (li >= lc) {
                        break;
                    }
                    u_elements = _mm512_set1_epi32(larray[li]);
                }
            }
        }
    } else {
        if (li + 3 < lc && ri + 3 < rc) {
            __m512i u_elements = _mm512_loadu_si512((__m512i *) (larray.data() + li));
            __m512i u_elements_per = _mm512_permutexvar_epi32(st, u_elements);
            __m512i v_elements = _mm512_loadu_si512((__m512i *) (rarray.data() + ri));
            // __m512i v_elements_per = _mm512_permutex_epi64(v_elements, 0b01000100);
            __m512i v_elements_per = _mm512_permutexvar_epi32(p4f128, v_elements);

            while (true) {
                __mmask16 mask = _mm512_cmpeq_epi32_mask(u_elements_per, v_elements_per);
                if (mask != 0x0000) {
                    // write back
                    _mm512_mask_compressstoreu_epi32(cur_back_ptr, mask, u_elements_per);
                    cur_back_ptr += _popcnt32(mask);
                }

                if (larray[li + 3] > rarray[ri + 3]) {
                    ri += 4;
                    if (ri + 3 >= rc) {
                        break;
                    }
                    v_elements = _mm512_loadu_si512((__m512i *) (rarray.data() + ri));
                    v_elements_per = _mm512_permutexvar_epi32(p4f128, v_elements);// _mm512_permutex_epi64(v_elements, 0b01000100);;
                } else if (larray[li + 3] < rarray[ri + 3]) {
                    li += 4;
                    if (li + 3 >= lc) {
                        break;
                    }
                    u_elements = _mm512_loadu_si512((__m512i *) (larray.data() + li));
                    u_elements_per = _mm512_permutexvar_epi32(st, u_elements);
                } else {
                    li += 4;
                    ri += 4;
                    if (li + 3 >= lc || ri + 3 >= rc) {
                        break;
                    }
                    u_elements = _mm512_loadu_si512((__m512i *) (larray.data() + li));
                    u_elements_per = _mm512_permutexvar_epi32(st, u_elements);
                    v_elements = _mm512_loadu_si512((__m512i *) (rarray.data() + ri));
                    v_elements_per = _mm512_permutexvar_epi32(p4f128, v_elements);//_mm512_permutex_epi64(v_elements, 0b01000100);
                }
            }
        }
    }

    cn.assign(cn.data(), cur_back_ptr);

    if (li < lc && ri < rc) {
        while (true) {
            while (larray[li] < rarray[ri]) {
                li += 1;
                if (li >= lc) {
                    return;
                }
            }
            while (larray[li] > rarray[ri]) {
                ri += 1;
                if (ri >= rc) {
                    return;
                }
            }
            if (larray[li] == rarray[ri]) {
                // write back
                //cn[cn_count++] = larray[li];
                cn.emplace_back(larray[li]);
                li += 1;
                ri += 1;
                if (li >= lc || ri >= rc) {
                    return;
                }
            }
        }
    }
    return;
}