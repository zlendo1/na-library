#ifndef NONLINEAREQUATIONSOLVER_H
#define NONLINEAREQUATIONSOLVER_H

#include <cmath>
#include <limits>
#include <stdexcept>
#include <functional>
#include <vector>
#include <complex>
#include <random>
#include <algorithm>

template <typename FunType>
bool BracketRoot(FunType f, double x0, double &a, double &b, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4) {
    if (hinit < 0 || hmax < 0 || lambda < 0) {
        throw std::domain_error("Invalid parameters");
    }

    a = x0;

    double f1 = f(a), f2, h = hinit;
    const double eps = std::numeric_limits<double>::epsilon();

    while (std::abs(h) < hmax) {
        b = a + h;

        f2 = f(b);

        while (!std::isfinite(f2)) {
            h /= 2 * (1 + lambda);

            if (std::abs(h) <= std::abs(a) * eps) {
                return false;
            }

            b = a + h;

            f2 = f(b);
        }

        if (f1 * f2 <= 0) {
            return true;
        }

        h *= lambda;

        a = b;

        f1 = f2;
    }

    a = x0;
    f1 = f(a);
    h = hinit;

    while (std::abs(h) < hmax) {
        b = a - h;

        f2 = f(b);

        while (!std::isfinite(f2)) {
            h /= 2 * (1 + lambda);

            if (std::abs(h) <= std::abs(a) * eps) {
                return false;
            }

            b = a - h;

            f2 = f(b);
        }

        if (f1 * f2 <= 0) {
            std::swap(a, b);

            return true;
        }

        h *= lambda;

        a = b;

        f1 = f2;
    }

    return false;
}

enum RegulaFalsiMode {Unmodified, Illinois, Slavic, IllinoisSlavic};

template <typename FunType>
double RegulaFalsiSolve(FunType f, double a, double b, RegulaFalsiMode mode = Slavic, double eps = 1e-10, int maxiter = 100) {
    if (f(a) * f(b) > 0) {
        throw std::range_error("Root must be bracketed");
    }

    if (eps < 0 || maxiter < 0) {
        throw std::domain_error("Invalid parameters");
    }

    std::function<double (double, FunType)> fi = [&] (double x, FunType f) { return f(x); };

    bool Illinois_Mode = mode == IllinoisSlavic || mode == Illinois, Slavic_Mode = mode == IllinoisSlavic || mode == Slavic;

    if (Slavic_Mode) {
        fi = [&] (double x, FunType) { return f(x) / (1 + std::abs(f(x))); };
    }

    int iter = 0;

    double f1 = fi(a, f), f2 = fi(b, f), c = a, cold = b;
    
    while (std::abs(c - cold) > eps && iter < maxiter) {
        cold = c;

        c = (a * f2 - b * f1) / (f2 - f1);

        double f3 = f(c);

        if (f3 == 0) {
            return c;
        }

        if (f1 * f3 < 0) {
            b = a;

            f2 = f1;
        } else if (Illinois_Mode) {
            f2 /= 2;
        }

        a = c;

        f1 = f3;

        ++iter;
    }

    if (iter == maxiter) {
        throw std::logic_error("Given accuracy has not achieved");
    }

    return c;
}

template<typename T>
int sgn(T number) {
    return number > 0 ? 1 : (number < 0 ? -1 : 0);
}

template <typename FunType>
double RiddersSolve(FunType f, double a, double b, double eps = 1e-10, int maxiter = 100) {
    if (f(a) * f(b) > 0) {
        throw std::range_error("Root must be bracketed");
    }

    if (eps < 0 || maxiter < 0) {
        throw std::domain_error("Invalid parameters");
    }
    
    int iter = 0;
    
    double f1 = f(a), f2 = f(b);

    while (std::abs(b - a) > eps && iter < maxiter) {
        double c = (a + b) / 2, f3 = f(c);

        if (f3 == 0) {
            return c;
        }

        double d = c + f3 * (c - a) * sgn(f1 - f2) / std::sqrt(f3 * f3 - f1 * f2), f4 = f(d);

        if (f4 == 0) {
            return d;
        }

        if (f3 * f4 <= 0) {
            a = c;
            b = d;

            f1 = f3;
            f2 = f4;
        } else if (f1 * f4 <= 0) {
            b = d;

            f2 = f4;
        } else {
            a = d;

            f1 = f4;
        }

        ++iter;
    }

    if (iter == maxiter) {
        throw std::logic_error("Given accuracy has not achieved");
    }

    return (a + b) / 2;
}

template <typename FunType1, typename FunType2>
double NewtonRaphsonSolve(FunType1 f, FunType2 fprim, double x0, double eps = 1e-10, double damping = 0, int maxiter = 100) {
    if (eps < 0 || maxiter < 0 || damping < 0 || damping >= 1) {
        throw std::domain_error("Invalid parameters");
    }
    
    int iter = 0;

    double delta;
    double v = f(x0), d = fprim(x0);

    do {
        if (std::abs(v) <= eps) {
            return x0;
        }

        delta = v / d;

        double w = v;

        v = f(x0 - delta);
        d = fprim(x0 - delta);

        while (damping && (!std::isfinite(v) || std::abs(v) > std::abs(w) || d == 0)) {
            delta *= damping;

            v = f(x0 - delta);
            d = fprim(x0 - delta);
        }

        x0 -= delta;
    } while (std::abs(x0) > eps && ++iter < maxiter);

    if (iter == maxiter || !std::isfinite(x0)) {
        throw std::logic_error("Convergence has not achieved");
    }

    return x0;
}

std::complex<double> RandomComplex(double min_r, double max_r, double min_i, double max_i) {
    std::random_device rd;
    std::default_random_engine rng(rd());

    #define normalize(number, min, max) std::fmod(double(number), (max - min)) + min

    return {normalize(rng(), min_r, max_r), normalize(rng(), min_i, max_i)};

    #undef normalize
}

template<typename T>
std::pair<std::complex<double>, bool> LaguerreRoot(std::vector<T> coefficients, int n, std::complex<double> x,  double eps, int maxtrials) {
    int trial = 0;

    std::complex<double> delta;

    do {
        std::complex<double> f = coefficients[n];
        std::complex<double> d = 0, s = 0;

        for (int i = n - 1; i >= 0; --i) {
            s = s * x + 2. * d;
            d = d * x + f;
            f = f * x + coefficients[i];
        }

        if (std::abs(f) <= eps) {
            return {x, true};
        }

        std::complex<double> r = std::sqrt(double(n - 1) * (double(n - 1) * d * d - double(n) * f * s));

        if (std::abs(d + r) > std::abs(d - r)) {
            delta = double(n) * f / (d + r);
        } else {
            delta = double(n) * f / (d - r);
        }

        x -= delta;
    } while (std::abs(delta) > eps && ++trial < maxtrials);

    if (std::abs(delta) <= eps) {
        return {x, true};
    }

    return {x, false};
}

std::vector<std::complex<double>> PolyRoots(std::vector<std::complex<double>> coefficients, double eps = 1e-10, int maxiters = 100, int maxtrials = 10) {
    if (eps < 0 || maxiters < 0 || maxtrials < 0) {
        throw std::domain_error("Invalid parameters");
    }

    int n = coefficients.size();

    std::vector<std::complex<double>> roots(n - 1);
    std::vector<std::complex<double>> original_coeff = coefficients;

    for (int i = n - 1; i >= 1; --i) {
        int iter = 0;
        bool c = false;

        std::complex<double> x;

        while (!c && iter < maxiters) {
            x = RandomComplex(-10, 10, -10, 10);

            std::tie(x, c) = LaguerreRoot(coefficients, i, x, eps, maxtrials);

            ++iter;
        }

        std::tie(x, c) = LaguerreRoot(original_coeff, n - 1, x, eps, maxtrials);

        if (!c) {
            throw std::logic_error("Convergence has not achieved");
        }

        if (std::abs(imag(x)) < eps) {
            x = real(x);
        }

        roots[i - 1] = x;

        std::complex<double> v = coefficients[i];

        for (int j = i - 1; j >= 0; --j) {
            std::complex<double> w = coefficients[j];

            coefficients[j] = v;

            v = w + x * v;
        }
    }

    std::sort(roots.begin(), roots.end(), [] (std::complex<double> a, std::complex<double> b) { return a.real() < b.real() || (a.real() == b.real() && a.imag() < b.imag()); });

    return roots;
}

std::vector<std::complex<double>> PolyRoots(std::vector<double> coefficients, double eps = 1e-10, int maxiters = 100, int maxtrials = 10) {
    if (eps < 0 || maxiters < 0 || maxtrials < 0) {
        throw std::domain_error("Invalid parameters");
    }

    int n = coefficients.size();

    std::vector<std::complex<double>> roots(n - 1);
    std::vector<double> original_coeff = coefficients;
    
    for (int i = n - 1; i >= 1;) {
        int iter = 0;
        bool c = false;

        std::complex<double> x;

        while (!c && iter++ < maxiters) {
            x = RandomComplex(-10, 10, -10, 10);

            std::tie(x, c) = LaguerreRoot(coefficients, i, x, eps, maxtrials);
        }

        std::tie(x, c) = LaguerreRoot(original_coeff, n - 1, x, eps, maxtrials);

        if (!c) {
            throw std::logic_error("Convergence has not achieved");
        }

        if (std::abs(imag(x)) <= eps) {
            x = real(x);

            roots[i - 1] = x;

            double v = coefficients[i];

            for (int j = i - 1; j >= 0; --j) {
                double w = coefficients[j];

                coefficients[j] = v;

                v = w + real(x) * v;
            }

            i -= 1;
        } else {
            roots[i - 1] = x;
            roots[i - 2] = std::conj(x);

            double a = 2 * real(x);
            double b = std::norm(x);

            double u = coefficients[i];
            double v = coefficients[i - 1] + a * u;

            for (int j = i - 2; j >= 0; --j) {
                double w = coefficients[j];

                coefficients[j] = u;

                u = v;
                v = w + a * v - b * coefficients[j];
            }

            i -= 2;
        }
    }

    std::sort(roots.begin(), roots.end(), [] (std::complex<double> a, std::complex<double> b) { return a.real() < b.real() || (a.real() == b.real() && a.imag() < b.imag()); });

    return roots;
}

#endif