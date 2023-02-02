#include <cmath>
#include <gtest/gtest.h>
#include <cosmology/utility/spline.h>
#include <cosmology/utility/math.h>

double calculateRelativeError(const std::vector<double> &base, const std::vector<double> &result) {
    double error = 0;

    for (size_t i = 0; i < base.size(); ++i) {
        double relative = std::abs(base[i] - result[i]) / (std::abs(base[i]) + std::abs(result[i]));
        if (error < relative) {
            error = relative;
        }
    }

    return error;
}

std::vector<double> map(const std::vector<double> &vector, const std::function<double(double)> &func) {
    std::vector<double> result;
    result.reserve(vector.size());

    for (const auto &element: vector) {
        result.push_back(func(element));
    }

    return result;
}

TEST(Spline, Cubic) {
    auto poly = [](double x) { return x * x * x - 3 * x * x + 5 * x + 8; };
    std::vector<std::pair<double, double>> points;

    for (double x: linSpace(-3.0, 3.0, 6)) {
        points.emplace_back(x, poly(x));
    }

    auto spline = SplineBuilder().buildFromPoints(std::move(points));
    auto test_points = linSpace(-3.0, 3.0, 100);
    ASSERT_LE(calculateRelativeError(map(test_points, poly),
                                     map(test_points, std::bind(&Spline::evaluate, spline, std::placeholders::_1))),
              0.07);
}

TEST(Spline, Exp) {
    auto func = [](double x) { return std::exp(x);};
    std::vector<std::pair<double, double>> points;

    for (double x: linSpace(-3.0, 3.0, 20)) {
        points.emplace_back(x, func(x));
    }

    auto spline = SplineBuilder().buildFromPoints(std::move(points));
    auto test_points = linSpace(-3.0, 3.0, 100);
    ASSERT_LE(calculateRelativeError(map(test_points, func),
                                     map(test_points, std::bind(&Spline::evaluate, spline, std::placeholders::_1))),
              0.003);
}

int main() {
    ::testing::InitGoogleTest();
    return RUN_ALL_TESTS();
}