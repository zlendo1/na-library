#ifndef LIMIT_H
#define LIMIT_H

#include <limits>
#include <utility>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <vector>

template <typename FunType>
std::pair<double, bool> Limit(FunType f, double x0, double h = 0, double eps = 1e-8, double nmax = 20) {
    if (eps <= 0 || nmax < 3 || nmax > 30) {
        throw std::domain_error("Invalid parameters");
    }
    
    double y_old = std::numeric_limits<double>::infinity();

    bool inf_flag = false, neg_inf_flag = false;
    
    if (std::isinf(x0)) {
        if (x0 < 0) {
            neg_inf_flag = true;
        }

        x0 = 0;

        inf_flag = true;
    }

    h = h != 0 ? h : 1e-3L * std::max(1., std::abs(x0));
    h = !neg_inf_flag ? h : -h;

    std::vector<double> y(nmax);

    for (int i = 0; i < nmax; ++i) {
        y[i] = !inf_flag ? f(x0 + h) : f(1 / (x0 + h));

        int p = 2;

        for (int k = i - 1; k >= 0; --k) {
            y[k] = ((long double) p * y[k + 1] - y[k]) / (p - 1);

            p *= 2;
        }

        if (std::abs(y.front() - y_old) < eps) {
            return {y.front(), true};
        }

        y_old = y.front();

        h /= 2L;
    }

    return {y.front(), false};
}

#endif