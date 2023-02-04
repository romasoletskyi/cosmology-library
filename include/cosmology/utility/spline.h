#pragma once

#include <vector>
#include <ostream>
#include <istream>

struct CubicPolynomial {
    double a, b, c, d;

    double evaluate(double t) const {
        return a + b * t + c * t * t + d * t * t * t;
    }
};

std::ostream &operator<<(std::ostream &stream, const CubicPolynomial &poly);

std::istream &operator>>(std::istream &stream, CubicPolynomial &poly);

class Spline {
public:
    Spline() {}

    Spline(const std::vector<double> &knots, const std::vector<CubicPolynomial> &polys) : knots_(knots),
                                                                                          polys_(polys) {}

    double evaluate(double t) const;

    const std::vector<double> &getKnots() const;

    const std::vector<CubicPolynomial> &getPolynomials() const;

private:
    friend class SplineBuilder;

    std::vector<double> knots_;
    std::vector<CubicPolynomial> polys_;
};

class SplineBuilder {
public:
    Spline buildFromPoints(std::vector<std::pair<double, double>> points);
};