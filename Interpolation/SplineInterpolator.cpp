#ifndef SPLINEINTERPOLATOR_H
#define SPLINEINTERPOLATOR_H

#include "AbstractInterpolator.cpp"

#define x(i) points[i].first
#define y(i) points[i].second

class SplineInterpolator : public AbstractInterpolator {
private:
    std::vector<double> r;
    std::vector<double> s;
    std::vector<double> q;
public:
    SplineInterpolator(const std::vector<std::pair<double, double>> &data);
    
    double operator()(double x) const override;
};

SplineInterpolator::SplineInterpolator(const std::vector<std::pair<double, double>> &data) :AbstractInterpolator(data), r(data.size()), s(data.size() - 1), q(data.size() - 1) {
    if (points.size() < 3) {
        throw std::domain_error("Invalid data set");
    }

    const int n = points.size();

    std::vector<long double> alpha(n - 2);

    r.front() = 0;
    r.back() = 0;

    for (int i = 1; i < n - 1; ++i) {
        alpha[i - 1] = 2 * ((long double) x(i + 1) - x(i - 1));

        r[i] = 3 * (((long double) y(i + 1) - y(i)) / ((long double) x(i + 1) - x(i)) - ((long double) y(i) - y(i - 1)) / ((long double) x(i) - x(i - 1)));
    }

    for (int i = 1; i < n - 2; ++i) {
        long double mi = ((long double) x(i + 1) - x(i)) / alpha[i - 1];

        alpha[i] -= mi * ((long double) x(i + 1) - x(i));

        r[i + 1] -= mi * r[i];
    }

    r[n - 2] /= alpha.back();

    for (int i = n - 3; i >= 1; --i) {
        r[i] = (r[i] - ((long double) x(i + 1) - x(i)) * r[i + 1]) / alpha[i - 1];
    }

    for (int i = 0; i < n - 1; ++i) {
        long double delta_x = (long double) x(i + 1) - x(i);

        s[i] = ((long double) r[i + 1] - r[i]) / (3 * delta_x);

        q[i] = ((long double) y(i + 1) - y(i)) / delta_x - delta_x * (r[i + 1] + 2L * r[i]) / 3;
    }
}
    
double SplineInterpolator::operator()(double x) const {
    int index_upper = Locate(x);

    if (index_upper <= 0) {
        ++index_upper;
    }

    if (index_upper >= points.size()) {
        --index_upper;
    }

    const int i = index_upper - 1;
    long double t = x - x(i);

    return y(i) + t * (q[i] + t * (r[i] + t * s[i]));
}

#undef x
#undef y

#endif