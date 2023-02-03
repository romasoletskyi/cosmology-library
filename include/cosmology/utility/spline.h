#include <vector>

struct CubicPolynomial {
    double a, b, c, d;

    double evaluate(double t) const {
        return a + b * t + c * t * t + d * t * t * t;
    }
};

class Spline {
public:
    double evaluate(double t) const;

    const std::vector<double>& getKnots() const;

private:
    friend class SplineBuilder;

    std::vector<double> knots_;
    std::vector<CubicPolynomial> polys_;
};

class SplineBuilder {
public:
    Spline buildFromPoints(std::vector<std::pair<double, double>> points);
};