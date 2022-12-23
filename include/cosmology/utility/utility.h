#pragma once

#include <vector>
#include <unordered_map>

template<class T>
std::unordered_map<T, int> createVectorToIndex(const std::vector<T> &vector) {
    std::unordered_map<T, int> map;
    for (int i = 0; i < vector.size(); ++i) {
        map[vector[i]] = i;
    }
    return map;
}
