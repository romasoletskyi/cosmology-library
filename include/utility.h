#include <vector>

template<class T>
std::vector<T> linSpace(T begin, T end, size_t pointsNum) {
    std::vector<T> vector;
    vector.reserve(pointsNum + 1);

    auto step = (end - begin) / pointsNum;
    auto value = begin;

    for (size_t i = 0; i <= pointsNum; ++i) {
        vector.push_back(value);
        value += step;
    }

    return vector;
}
