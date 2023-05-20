#ifndef DIFFERENTIALSOLVER_H
#define DIFFERENTIALSOLVER_H

#include <stdexcept>
#include <vector>
#include <cmath>

template <typename FunType>
double RK4Step(FunType f, double x, double y, double h) {
    double K1 = f(x, y);
    double K2 = f(x + h / 2, y + h * K1 / 2);
    double K3 = f(x + h / 2, y + h * K2 / 2);
    double K4 = f(x + h, y + h * K3);

    return y + h * (K1 + 2 * K2 + 2 * K3 + K4) / 6;
}

template <typename FunType>
std::vector<std::pair<double, double>> RK4Integrator(FunType f, double x0, double y0, double xmax, double h, double eps = 1e-8, bool adaptive = false) {
    if (h == 0 || eps < 0) {
        throw std::domain_error("Invalid paramaters");
    }

    std::vector<std::pair<double, double>> solutions;

    double x = x0;
    double y = y0;

    while (h > 0 && x <= xmax + eps || h < 0 && x >= xmax - eps) {
        if (adaptive) {
            double u = RK4Step(f, x, y, h / 2);
            double v = RK4Step(f, x + h / 2, u, h / 2);
            double w = RK4Step(f, x, y, h);
            double delta = std::abs(w - v) / std::abs(h);

            if (delta <= eps) {
                if (h > 0 && x + h > xmax || h < 0 && x + h < xmax) {
                    h = xmax - x;

                    continue;
                } else {
                    x += h;
                }

                y = v;

                solutions.emplace_back(x, y);
            }

            h *= std::min(5., 0.9 * std::pow(eps / delta, 0.25));

        } else {
            solutions.emplace_back(x, y);

            y = RK4Step(f, x, y, h);

            x += h;
        }
    }

    return solutions;
}

template <typename FunType>
std::vector<std::pair<double, std::vector<double>>> RK4SystemIntegrator(FunType f, double x0, std::vector<double> y0, double xmax, double h) {
    if (h == 0) {
        throw std::domain_error("Invalid paramaters");
    }

    std::vector<std::pair<double, std::vector<double>>> solutions;

    const int n = y0.size();

    double x = x0;
    std::vector<double> y = y0;

    solutions.emplace_back(x, y);

    while (h > 0 && x <= xmax || h < 0 && x >= xmax) {
        std::vector<std::vector<double>> K(4, std::vector<double>(n));
        std::vector<double> u(n), v(n);

        K[0] = f(x, y);

        for (int k = 0; k < n; ++k) {
            u[k] = y[k] + h * K[0][k] / 2;
        }

        x += h / 2;

        K[1] = f(x, u);

        for (int k = 0; k < n; ++k) {
            v[k] = y[k] + h * K[1][k] / 2;
        }

        K[2] = f(x, v);

        for (int k = 0; k < n; ++k) {
            u[k] = y[k] + h * K[2][k];
        }

        x += h / 2;

        K[3] = f(x, u);

        for (int k = 0; k < n; ++k) {
            y[k] += h * (K[0][k] + 2 * K[1][k] + 2 * K[2][k] + K[3][k]) / 6;
        }

        solutions.emplace_back(x, y);
    }

    return solutions;
}

#endif