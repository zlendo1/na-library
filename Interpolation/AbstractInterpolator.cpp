#ifndef ABSTRACTINTERPOLATOR_H
#define ABSTRACTINTERPOLATOR_H

#include <iterator>
#include <stdexcept>
#include <vector>
#include <utility>
#include <algorithm>
#include <functional>

#define x(i) points[i].first

class AbstractInterpolator {
private:
protected:
    using Point = std::pair<double, double>;
    using VectorPoint = std::vector<Point>;

    VectorPoint points;

public:
    int Locate(double x) const;

private:
    mutable int cached_index = 0;

public:
    AbstractInterpolator(const std::vector<std::pair<double, double>> &data);

    virtual double operator()(double x) const = 0;
};

AbstractInterpolator::AbstractInterpolator(const std::vector<std::pair<double, double>> &data) {
    if (data.size() < 2) {
        throw std::domain_error("Invalid data set");
    }

    points = data;

    std::sort(points.begin(), points.end(), [] (const Point &a, const Point &b) { return a.first < b.first; });

    for (int i = 0; i < data.size() - 1; ++i) {
        if (x(i) == x(i + 1)) {
            throw std::domain_error("Invalid data set");
        }
    }
}

int AbstractInterpolator::Locate(double x) const {
    if (cached_index && points[cached_index - 1].first > x && points[cached_index].first <= x) {
        return cached_index;
    }

    if (x <= points.front().first) {
        return 0;
    }

    if (x > points.back().first) {
        return points.size();
    }

    auto iterator = std::upper_bound(points.begin(), points.end(), x, [] (double a, const Point &b) { return a <= b.first; });

    cached_index =std::distance(points.begin(), iterator);

    return cached_index;
}

#undef x

#endif