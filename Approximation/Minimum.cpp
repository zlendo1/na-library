#ifndef MINIMUM_H
#define MINIMUM_H

#include <stdexcept>
#include <limits>
#include <cmath>

template <typename FunType>
bool BracketMinimum(FunType f, double x0, double &a, double &b, double eps, double hinit, double hmax, double lambda) {
    if (hinit < 0 || hmax < 0 || lambda < 0) {
        throw std::domain_error("Invalid parameters");
    }

    a = x0;

    double f1 = f(a), f2, h = hinit, c, f3;

    while (std::abs(h) < hmax) {
        b = a + h;
        c = a + h / 2;

        f2 = f(b);
        f3 = f(c);

        while (!std::isfinite(f2)) {
            h /= 2 * (1 + lambda);

            if (std::abs(h) <= std::abs(a) * eps) {
                return false;
            }

            b = a + h;
            c = a + h / 2;

            f2 = f(b);
            f3 = f(c);
        }

        if (f1 >= f3 && f2 >= f3) {
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
        c = a - h / 2;

        f2 = f(b);
        f3 = f(c);

        while (!std::isfinite(f2)) {
            h /= 2 * (1 + lambda);

            if (std::abs(h) <= std::abs(a) * eps) {
                return false;
            }

            b = a - h;
            c = a - h / 2;

            f2 = f(b);
            f3 = f(c);
        }

        if (f1 >= f3 && f2 >= f3) {
            return true;
        }

        h *= lambda;

        a = b;

        f1 = f2;
    }

    return false;
}

template <typename FunType>
double FindMinimum(FunType f, double x0, double eps = 1e-8, double hinit = 1e-5, double hmax = 1e10, double lambda = 1.4) {
    if (eps < 0 || hinit < 0 || hmax < 0 || lambda < 0) {
        throw std::domain_error("Invalid parameters");
    }

    double a, b;

    if (!BracketMinimum(f, x0, a, b, eps, hinit, hmax, lambda)) {
        throw std::logic_error("Minimum has not found");
    }

    const double fi = (1 + sqrt(5)) / 2;
    double c = b - (b - a) / fi;
    double d = a + (b - a) / fi;
    double u = f(c);
    double v = f(d);

    while (std::abs(b - a) > eps) {
        if (u < v) {
            b = d;
            d = c;

            c = a + (c - a) / fi;

            v = u;
            u = f(c);
        } else {
            a = c;
            c = d;

            d = b - (b - d) / fi;

            u = v;
            v = f(d);
        }
    }

    return (a + b) / 2;
}

#endif