#include <random>

#include <gtest/gtest.h>
#include <cosmology/utility/math.h>

std::vector<double> map(const std::vector<double> &grid, const std::function<double(double)> &function) {
    std::vector<double> result;
    for (double x: grid) {
        result.push_back(function(x));
    }
    return result;
}

void checkDerivativeOnParabola(const std::vector<double>& time) {
    auto array = map(time, [](double t) { return t * t; });

    auto result = map(time, [](double t) { return 2 * t; });
    double test[result.size()];
    differentiate(time.data(), array.data(), test, result.size());

    for (size_t i = 0; i < result.size(); ++i) {
        ASSERT_NEAR(test[i], result[i], 1e-9);
    }
}

TEST(Derivative, Basic) {
    checkDerivativeOnParabola(linSpace<double>(-0.5, 1.0, 17));
}

TEST(Derivative, VariableInterval) {
    std::mt19937 gen;
    std::uniform_real_distribution<> uni(0.05, 0.15);

    std::vector<double> grid;
    double time = -0.7;

    while (time < 2) {
        grid.push_back(time);
        time += uni(gen);
    }

    checkDerivativeOnParabola(grid);
}

int main() {
    ::testing::InitGoogleTest();
    //::testing::GTEST_FLAG(filter) = "Derivative*";
    return RUN_ALL_TESTS();
}
