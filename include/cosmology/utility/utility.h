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

template<class T>
std::vector<std::pair<T, T>> zip(const std::vector<T>& lhs, const std::vector<T>& rhs) {
    size_t length = std::min(lhs.size(), rhs.size());
    std::vector<T> result;
    result.reserve(length);

    for (size_t i = 0; i < length; ++i) {
        result.emplace_back(lhs[i], rhs[i]);
    }

    return result;
}
