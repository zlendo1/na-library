#ifndef LINEARINTERPOLATOR_H
#define LINEARINTERPOLATOR_H

#include "AbstractInterpolator.cpp"

class LinearInterpolator : public AbstractInterpolator {
public:
    LinearInterpolator(const std::vector<std::pair<double, double>> &data);

    double operator()(double x) const override;
};

LinearInterpolator::LinearInterpolator(const std::vector<std::pair<double, double>> &data) : AbstractInterpolator(data) {}
    
double LinearInterpolator::operator()(double x) const {
    int index_upper = Locate(x);

    if (index_upper <= 0) {
        ++index_upper;
    }

    if (index_upper >= points.size()) {
        --index_upper;
    }

    const Point &lower_limit = points[index_upper - 1];
    const Point &upper_limit = points[index_upper];

    if (upper_limit.first == x) {
        return upper_limit.second;
    }

    long double range_size = (long double) upper_limit.first - lower_limit.first;
    long double x_long = x;

    return (upper_limit.first - x_long) / range_size * lower_limit.second +
           (x_long - lower_limit.first) / range_size * upper_limit.second;
}

#endif