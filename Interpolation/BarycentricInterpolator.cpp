#ifndef BARYCENTRICINTERPOLATOR_H
#define BARYCENTRICINTERPOLATOR_H

#include "AbstractInterpolator.cpp"
#include <algorithm>
#include <stdexcept>

#define x(i) points[i].first
#define y(i) points[i].second

class BarycentricInterpolator : public AbstractInterpolator {
private:
    std::vector<double> w;

public:
    BarycentricInterpolator(const std::vector<std::pair<double, double>> &data, int order);

    double operator()(double x) const override;

    std::vector<double> GetWeights() const;
};

BarycentricInterpolator::BarycentricInterpolator(const std::vector<std::pair<double, double>> &data, int order) : AbstractInterpolator(data), w(data.size(), 0) {
    const int n = points.size();
    
    if (order < 0 || order > n) {
        throw std::domain_error("Invalid order");
    }

    const int d = order;

    for (int i = 1; i <= n; ++i) {
        long double p = 1;

        for (int k = std::max(1, i - d); k <= std::min(i, n - d); ++k) {
            for (int j = k; j <= k + d; ++j) {
                if (i != j) {
                    p /= (long double) x(i - 1) - x(j - 1);
                }
            }

            if (k % 2 == 0) {
                p = -p;
            }
        }

        w[i - 1] += p;
    }
}

double BarycentricInterpolator::operator()(double x) const {
    long double p = 0, q = 0;
    const int n = points.size();

    for (int i = 0; i < n; ++i) {
        if (x == x(i)) {
            return y(i);
        }

        long double u = w[i] / (x - x(i));

        p += u * y(i);
        q += u;
    }

    return p / q;
}

std::vector<double> BarycentricInterpolator::GetWeights() const {
    return w;
}

#undef x
#undef y

#endif