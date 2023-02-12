#pragma once

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

template<class Time>
std::tuple<Time, Time, Time> calculateDifferentiationCoefficients(Time left, Time right) {
    return {-right / left / (left + right), (right - left) / left / right, left / right / (left + right)};
}

template<class Time, class Value>
void differentiate(const Time *time, const Value *array, Value *result, size_t size) {
    {
        auto [alpha, beta, gamma] = calculateDifferentiationCoefficients(time[0] - time[1], time[2] - time[0]);
        result[0] = alpha * array[1] + beta * array[0] + gamma * array[2];
    }

    for (size_t i = 1; i < size - 1; ++i) {
        auto [alpha, beta, gamma] = calculateDifferentiationCoefficients(time[i] - time[i - 1], time[i + 1] - time[i]);
        result[i] = alpha * array[i - 1] + beta * array[i] + gamma * array[i + 1];
    }

    {
        auto [alpha, beta, gamma] = calculateDifferentiationCoefficients(time[size - 1] - time[size - 2],
                                                                         time[size - 3] - time[size - 1]);
        result[size - 1] = alpha * array[size - 2] + beta * array[size - 1] + gamma * array[size - 3];
    }
}

template<class Time, class Value>
Value integrate(const Time *time, const Value *array, size_t size) {
    Value sum = 0;
    for (int i = 0; i < size - 1; ++i) {
        sum += (time[i + 1] - time[i]) * (array[i + 1] + array[i]);
    }
    return 0.5 * sum;
}
