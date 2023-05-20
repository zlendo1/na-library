#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <stdexcept>
#include <utility>
#include <vector>
#include <list>
#include <cmath>

void CheckCommonParamaters(double eps, int nmax, int nmin) {
    if (eps < 0 || nmin < 0 || nmax < 0 || nmax < nmin) {
        throw std::domain_error("Bad parameter");
    }
}

template <typename FunType>
std::pair<double, bool> RombergIntegration(FunType f, double a, double b, double eps = 1e-8, int nmax = 1000000, int nmin = 50) {
    CheckCommonParamaters(eps, nmax, nmin);

    int N = 2;

    double h = (b - a) / N;
    double s = (f(a) + f(b)) / 2;
    double Iold = s;

    std::list<double> I;

    while (N <= nmax) {
        for (int j = 1; j <= N / 2; ++j) {
            s += f(a + (2 * j - 1) * h);
        }

        I.push_back(h * s);

        double Iprev = I.back();
        
        unsigned long long int p = 4;

        for (auto iter = ++I.rbegin(); iter != I.rend(); ++iter) {
            *iter = (p * Iprev - *iter) / (p - 1);

            Iprev = *iter;

            p *= 4;
        }

        if (N >= nmin && std::abs(I.front() - Iold) <= eps) {
            return std::pair<double, bool>(I.front(), true);
        }

        Iold = I.front();

        h /= 2;
        N *= 2;
    }

    return std::pair<double, bool>(I.front(), false);
}

template <typename FunType>
std::pair<double, bool> TanhSinhIntegration(FunType f, double a, double b, double eps = 1e-8, int nmax = 1000000, int nmin = 20, double range = 3.5) {
    CheckCommonParamaters(eps, nmax, nmin);

    if (range < 0) {
        throw std::domain_error("Bad parameter");
    }

    int N = 2;

    const double PI = 4 * std::atan(1L);

    double h = 2 * range / N;
    double p = (b + a) / 2;
    double q = (b - a) / 2;
    double s = 0;
    double Iold = 0;
    double I = 0;

    while (N < nmax) {
        for (int i = 1; i <= N / 2; ++i) {
            double t = (2 * i - 1) * h - range;
            double u = PI * std::sinh(t) / 2;
            double v = f(p + q * std::tanh(u));

            if (std::isfinite(v)) {
                s += q * PI * std::cosh(t) * v / (2 * std::cosh(u) * std::cosh(u));
            }
        }

        I = h * s;

        if (N >= nmin && std::abs(I - Iold) <= eps) {
            return std::pair<double, bool>(I, true);
        }

        Iold = I;

        N *= 2;
        h /= 2;
    }

    return std::pair<double, bool>(I, false);
}

template <typename FunType>
std::pair<double, bool> AdaptiveAux(FunType f, double a, double b, double eps, double f1, double f2, double f3, int depth) {
    f1 = std::isfinite(f1) ? f1 : 0;
    f2 = std::isfinite(f2) ? f2 : 0;
    
    double c = (a + b) / 2;

    if (!std::isfinite(f1) || !std::isfinite(f2)) {
        return AdaptiveAux(f, a, b, eps, f(a), f(b), f(c), depth);
    }

    double I1 = (b - a) * (f1 + 4 * f3 + f2) / 6;
    
    double f4 = f((a + c) / 2);
    double f5 = f((c + b) / 2);

    double I2 = (b - a) * (f1 + 4 * f4 + 2 * f3 + 4 * f5 + f2) / 12;

    if (std::abs(I1 - I2) <= eps) {
        return std::pair<double, bool>(I2, true);
    }

    if (depth <= 0) {
        return std::pair<double, bool>(I2, false);
    }

    auto ax1 = AdaptiveAux(f, a, c, eps, f1, f3, f4, depth - 1);
    auto ax2 = AdaptiveAux(f, c, b, eps, f3, f2, f5, depth - 1);

    return std::pair<double, int>(ax1.first + ax2.first, ax1.second && ax2.second);
}

template <typename FunType>
std::pair<double, bool> AdaptiveIntegration(FunType f, double a, double b, double eps = 1e-10, int maxdepth = 30, int nmin = 1) {
    if (eps < 0 || maxdepth < 0 || nmin < 0) {
        throw std::domain_error("Bad parameter");
    }
    
    double s = 0;
    double h = (b - a) / nmin;

    bool isValid = true;

    for (int i = 1; i <= nmin; ++i) {
        auto currentPair = AdaptiveAux(f, a, a + h, eps, f(a), f(a + h), f(a + h / 2), maxdepth);
        
        s += currentPair.first;
        isValid = !currentPair.second ? false : isValid;

        a += h;
    }

    return std::pair<double, bool>(s, isValid);
}

#endif