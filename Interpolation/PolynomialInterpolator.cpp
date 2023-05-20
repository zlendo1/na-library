#ifndef POLYNOMIALINTERPOLATOR_H
#define POLYNOMIALINTERPOLATOR_H

#include "AbstractInterpolator.cpp"
#include <algorithm>
#include <stdexcept>
#include <vector>

#define x(i) points[i].first
#define y(i) points[i].second

class PolynomialInterpolator : public AbstractInterpolator {
private:
    std::vector<double> newton_coefficients;
    std::vector<double> last_diagonal;

    void ComputeCoefficients(int j_bound);

public:
    PolynomialInterpolator(const std::vector<std::pair<double, double>> &data);

    double operator()(double x) const override;

    void AddPoint(const std::pair<double, double> &p);

    std::vector<double> GetCoefficients() const;
};

PolynomialInterpolator::PolynomialInterpolator(const std::vector<std::pair<double, double>> &data) : AbstractInterpolator(data), newton_coefficients(data.size()), last_diagonal(data.size()) {
    const int n = newton_coefficients.size();
    std::vector<double> &y = newton_coefficients;
    
    for (int i = 0; i < n; ++i) {
        y[i] = y(i);
    }

    last_diagonal.front() = y.back();

    for (int j = 0; j < n - 1; ++j) {
        for (int i = n - 1; i > j; --i) {
            y[i] = ((long double) y[i] - y[i - 1]) / ((long double) x(i) - x(i - j - 1));
        }

        last_diagonal[j + 1] = y.back();
    }
}

double PolynomialInterpolator::operator()(double x) const {
    long double f_x = newton_coefficients.back();

    const std::vector<double> &q = newton_coefficients;
    const int n = newton_coefficients.size();

    for (int i = n - 2; i >= 0; --i) {
        f_x = f_x * ((long double) x - x(i)) + q[i];
    }

    return f_x;
}

void PolynomialInterpolator::AddPoint(const std::pair<double, double> &p) {
    if (std::any_of(points.begin(), points.end(), [&] (const Point &x) { return p.first == x.first; })) {
        throw std::domain_error("Invalid point");
    }

    points.push_back(p);

    const int n = newton_coefficients.size();
    std::vector<double> &y = newton_coefficients;

    y.push_back(p.second);

    for (int j = 0; j < n; ++j) {
        double old_y = y[n];

        y[n] = ((long double) y[n] - last_diagonal[j]) / ((long double) x(n) - x(n - j - 1));

        last_diagonal[j] = old_y; 
    }

    last_diagonal.push_back(y.back());
}

std::vector<double> PolynomialInterpolator::GetCoefficients() const {
    const int n = points.size();
    
    std::vector<double> p(n, 0), w(n + 1, 1);

    for (int i = 1; i <= n; ++i) {
        w[i] = w[i - 1];

        for (int j = i - 1; j > 0; --j) {
            w[j] = w[j - 1] - x(i - 1) * w[j];
        }

        w.front() *= - x(i - 1);
    }

    for (int i = 0; i < n; ++i) {
        long double alpha = 1;

        for (int j = 0; j < n; ++j) {
            if (i != j) {
                alpha *= ((long double) x(i) - x(j));
            }
        }

        alpha = y(i) / alpha;

        std::vector<double> v(w);

        for (int j = n - 1; j >= 0; --j) {
            v[j] += (long double) x(i) * v[j + 1];

            p[j] += alpha * v[j + 1];
        }
    }

    return p;
}

#undef x
#undef y

#endif