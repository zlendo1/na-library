#ifndef TRIGONOMETRICINTERPOLATOR_H
#define TRIGONOMETRICINTERPOLATOR_H

#include "AbstractInterpolator.cpp"
#include <cmath>
#include <limits>
#include <stdexcept>

#define x(i) points[i].first
#define y(i) points[i].second

class TrigonometricInterpolator : public AbstractInterpolator {
private:
    double OMEGA_HALF;

    double EvenInterpolation(double x) const;
    double OddInterpolation(double x) const;

public:
    TrigonometricInterpolator(const std::vector<std::pair<double, double>> &data);

    double operator()(double x) const override;
};

TrigonometricInterpolator::TrigonometricInterpolator(const std::vector<std::pair<double, double>> &data) : AbstractInterpolator(data) {
    long double y_f = points.front().second;
    long double y_b = points.back().second;
    
    if (std::abs(y_f - y_b) < std::numeric_limits<double>::epsilon() * std::abs(y_f + y_b)) {
        throw std::domain_error("Function is not periodic");
    }

    long double x_f = points.front().first;
    long double x_b = points.back().second;

    OMEGA_HALF = 4 * std::atan(1L) / (x_b - x_f);
}

double TrigonometricInterpolator::operator()(double x) const {
    const int n = points.size();
    
    if (n % 2 == 0) {
        return EvenInterpolation(x);
    } else {
        return OddInterpolation(x);
    }
}

double TrigonometricInterpolator::EvenInterpolation(double x) const {
    const int n = points.size();
    long double f = 0;

    for (int i = 0; i < n - 1; ++i) {
        long double mul = 1;

        for (int j = 0; j < n - 1; ++j) {
            if (i == j) {
                continue;
            }

            mul *= std::sin(OMEGA_HALF * ((long double) x - x(j))) / std::sin(OMEGA_HALF * ((long double) x(i) - x(j)));
        }

        f += y(i) * mul;
    }

    return f;
}

double TrigonometricInterpolator::OddInterpolation(double x) const {
    const int n = points.size();
    long double f = 0;

    for (int i = 0; i < n - 1; ++i) {
        long double mul = 1;
        long double alpha = 0;

        for (int j = 0; j < n - 1; ++j) {
            if (i == j) {
                continue;
            }

            mul *= std::sin(OMEGA_HALF * ((long double) x - x(j))) / std::sin(OMEGA_HALF * ((long double) x(i) - x(j)));

            alpha -= x(j);
        }

        long double odd_factor = std::sin(OMEGA_HALF * (x - alpha)) / std::sin(OMEGA_HALF * (x(i) - alpha));

        f += y(i) * odd_factor * mul;
    }

    return f;
}

#undef x
#undef y

#endif