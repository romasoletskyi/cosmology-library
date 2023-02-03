#include <algorithm>

#include <cosmology/utility/spline.h>

double Spline::evaluate(double t) const {
    size_t index = std::distance(knots_.begin(), std::lower_bound(knots_.begin(), knots_.end(), t));
    if (index > 0) {
        --index;
    }
    if (index == polys_.size()) {
        --index;
    }
    return polys_[index].evaluate(t - knots_[index]);
}

const std::vector<double> &Spline::getKnots() const {
    return knots_;
}

// algorithm taken from https://en.wikipedia.org/wiki/Spline_(mathematics)
Spline SplineBuilder::buildFromPoints(std::vector<std::pair<double, double>> points) {
    std::sort(points.begin(), points.end(), [](const auto &lhs, const auto &rhs) {
        return lhs.first < rhs.first;
    });

    Spline spline;
    spline.knots_.reserve(points.size());
    for (const auto &point: points) {
        spline.knots_.push_back(point.first);
    }

    size_t size = points.size() - 1;
    double a[size + 1];
    for (size_t i = 0; i < size + 1; ++i) {
        a[i] = points[i].second;
    }

    double b[size], d[size], h[size];
    for (size_t i = 0; i < size; ++i) {
        h[i] = points[i + 1].first - points[i].first;
    }

    double alpha[size];
    for (size_t i = 1; i < size; ++i) {
        alpha[i] = 3 * (a[i + 1] - a[i]) / h[i] - 3 * (a[i] - a[i - 1]) / h[i - 1];
    }

    double c[size + 1], l[size + 1], mu[size + 1], z[size + 1];

    l[0] = 1;
    mu[0] = z[0] = 0;

    for (size_t i = 1; i < size; ++i) {
        l[i] = 2 * (points[i + 1].first - points[i - 1].first) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    l[size] = 1;
    z[size] = c[size] = 0;

    for (int i = static_cast<int>(size) - 1; i >= 0; --i) {
        c[i] = z[i] - mu[i] * c[i + 1];
        b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3;
        d[i] = (c[i + 1] - c[i]) / 3 / h[i];
    }

    spline.polys_.reserve(size);
    for (size_t i = 0; i < size; ++i) {
        spline.polys_.emplace_back(CubicPolynomial{a[i], b[i], c[i], d[i]});
    }

    return spline;
}
