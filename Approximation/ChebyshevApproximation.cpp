#ifndef CHEBYSHEVAPPROXIMATION_H
#define CHEBYSHEVAPPROXIMATION_H

#include <cmath>
#include <vector>
#include <stdexcept>

class ChebyshevApproximation {
    std::vector<double> chebyshevCoefficients;

    int m;
    double xmin, xmax;

    ChebyshevApproximation(const std::vector<double> &v, double xmin, double xmax) : chebyshevCoefficients(v), m(v.size() - 1), xmin(xmin), xmax(xmax) {}

    template <typename FunType>
    void CalculateCoefficients(FunType f);

    bool isInRange(double x) const;

public:
    template <typename FunType>
    ChebyshevApproximation(FunType f, double xmin, double xmax, int n);

    void set_m(int m);
    void trunc(double eps);

    double operator()(double x) const;

    double derivative(double x) const;
    ChebyshevApproximation derivative() const;

    ChebyshevApproximation antiderivative() const;
    double integrate(double a, double b) const;
    double integrate() const;
};

template <typename FunType>
void ChebyshevApproximation::CalculateCoefficients(FunType f) {
    std::vector<double> &c = chebyshevCoefficients, w(m + 2), v(m + 1);

    long double PI = 4 * std::atan(1L);
    const int n = m;

    for (int i = 0; i < n + 2; ++i) {
        w[i] = std::cos(PI * i / (2 * n + 2));
    }

    for (int i = 0; i <= n / 2; ++i) {
        v[i] = f((xmin + xmax + (xmax - xmin) * w[2 * i + 1]) / 2);
    }

    for (int i = n / 2 + 1; i < n + 1; ++i) {
        v[i] = f((xmin + xmax - (xmax - xmin) * w[2 * n + 1 - 2 * i]) / 2);
    }

    for (int k = 0; k < n + 1; ++k) {
        long double s = 0;

        for (int i = 0; i < n + 1; ++i) {
            int p = (k * (2 * i + 1)) % (4 * n + 4);

            if (p > 2 * n + 2) {
                p = 4 * n + 4 - p;
            }

            if (p > n + 1) {
                s -= v[i] * w[2 * n + 2 - p];
            } else {
                s += v[i] * w[p];
            }
        }

        c[k] = 2 * s / (n + 1);
    }
}

template <typename FunType>
ChebyshevApproximation::ChebyshevApproximation(FunType f, double xmin, double xmax, int n) : m(n), chebyshevCoefficients(n + 1), xmin(xmin), xmax(xmax) {
    if (n < 1 || xmin >= xmax) {
        throw std::domain_error("Bad parameters");
    }

    CalculateCoefficients(f);
}

void ChebyshevApproximation::set_m(int m) {
    if (m <= 1 || m >= chebyshevCoefficients.size()) {
        throw std::domain_error("Bad order");
    }

    ChebyshevApproximation::m = m;
}

void ChebyshevApproximation::trunc(double eps) {
    if (eps < 0) {
        throw std::domain_error("Bad tolerance");
    }

    std::vector<double> &c = chebyshevCoefficients;

    int new_m = c.size() - 1;

    for (auto iter = c.crbegin(); iter < c.crend(); iter++) {
        if (eps <= std::abs(*iter)) {
            break;
        }

        --new_m;
    }

    if (new_m < 1) {
        throw std::domain_error("Bad tolerance");
    }

    m = new_m;
}

double ChebyshevApproximation::operator()(double x) const {
    if (!isInRange(x)) {
        throw std::domain_error("Bad argument");
    }

    const std::vector<double> &c = chebyshevCoefficients;
    
    double t = (2 * x - xmin - xmax) / (xmax - xmin);
    double p = 0;
    double q = c[m];

    for (int k = m - 1; k > 0; --k) {
        p = 2 * t * q - p + c[k];
        
        std::swap(p, q);
    }

    return t * q - p + c.front() / 2;
}

double ChebyshevApproximation::derivative(double x) const {
    if (!isInRange(x)) {
        throw std::domain_error("Bad argument");
    }

    const std::vector<double> &c = chebyshevCoefficients;
    
    double t = (2 * x - xmin - xmax) / (xmax - xmin);
    double p = 0;
    double q = c[m];

    for (int k = m - 1; k > 1; --k) {
        p = (2 * (k + 1) * t * q - (k + 2) * p) / k + c[k];
        
        std::swap(p, q);
    }

    return 2 * (c[1] - 3 * p + 4 * t * q) / (xmax - xmin);
}

ChebyshevApproximation ChebyshevApproximation::derivative() const {
    const double mi = 4 / (xmax - xmin);

    const std::vector<double> &c = chebyshevCoefficients;
    std::vector<double> cp(m);

    cp[m - 1] = mi * m * c[m];
    cp[m - 2] = mi * (m - 1) * c[m - 1];

    for (int k = m - 3; k >= 0; --k) {
        cp[k] = cp[k + 2] + mi * (k + 1) * c[k + 1];
    }

    return ChebyshevApproximation(cp, xmin, xmax);
}

ChebyshevApproximation ChebyshevApproximation::antiderivative() const {
    const double mi = (xmax - xmin) / 4;

    const std::vector<double> &c = chebyshevCoefficients;
    std::vector<double> ci(m + 2);

    ci[0] = 0;
    ci[m + 1] = mi * c[m] / (m + 1);
    ci[m] = mi * c[m - 1] / m;

    for (int k = m - 1; k > 0; --k) {
        ci[k] = mi * (c[k - 1] - c[k + 1]) / k;
    }

    return ChebyshevApproximation(ci, xmin, xmax);
}

double ChebyshevApproximation::integrate(double a, double b) const {
    if (!isInRange(a) || !isInRange(b)) {
        throw std::domain_error("Bad interval");
    }

    ChebyshevApproximation integral = antiderivative(); 

    return integral(b) - integral(a);
}

double ChebyshevApproximation::integrate() const {
    const std::vector<double> &c = chebyshevCoefficients;

    double sum = c[0] / 2;

    for (int k = 2 ; k <= m; k += 2) {
        sum += c[k] / (1 - k * k);
    }

    return (xmax - xmin) * sum;
}

bool ChebyshevApproximation::isInRange(double x) const {
    return xmin <= x && x <= xmax;
}

#endif