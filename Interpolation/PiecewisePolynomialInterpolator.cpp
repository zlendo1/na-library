#ifndef PIECEWISEPOLYNOMIALINTERPOLATOR_H
#define PIECEWISEPOLYNOMIALINTERPOLATOR_H

#include "AbstractInterpolator.cpp"
#include <stdexcept>

class PiecewisePolynomialInterpolator : public AbstractInterpolator {
private:
    int order;

public:
    PiecewisePolynomialInterpolator(const std::vector<std::pair<double, double>> &data, int order);
    
    double operator()(double x) const override;
};

PiecewisePolynomialInterpolator::PiecewisePolynomialInterpolator(const std::vector<std::pair<double, double>> &data, int order) : AbstractInterpolator(data) {
    if (order < 1 || order > points.size()) {
        throw std::domain_error("Invalid order");
    }

    PiecewisePolynomialInterpolator::order = order;
}
    
double PiecewisePolynomialInterpolator::operator()(double x) const {
    int index_upper = Locate(x);

    if (index_upper <= 0) {
        ++index_upper;
    }

    if (index_upper >= points.size()) {
        --index_upper;
    }

    auto lower_iterator = points.begin() + index_upper - 1;
    auto upper_iterator = points.begin() + index_upper;

    int k = order;

    while (true) {
        if (upper_iterator < points.end()) {
            ++upper_iterator;

            --k;
        }

        if (!k) {
            break;
        }

        if (lower_iterator > points.begin()) {
            --lower_iterator;

            --k;
        }

        if (!k) {
            break;
        }
    }

    VectorPoint new_points(lower_iterator, upper_iterator);

    const int n = new_points.size();
    std::vector<long double> y(n);

    #define x(i) new_points[i].first

    for (int i = 0; i < n; ++i) {
        y[i] = new_points[i].second;
    }

    long double long_x = x;

    for (int j = 0; j < n - 1; ++j) {
        for (int i = n - 1; i > j; --i) {
            y[i] = (y[i] * (long_x - x(i - j - 1)) - y[i - 1] * (long_x - x(i))) / ((long double) x(i) - x(i - j - 1));
        }
    }

    #undef x

    return y.back();
}

#endif