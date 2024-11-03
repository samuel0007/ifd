#ifndef IFD_UTILS_HPP
#define IFD_UTILS_HPP

#include <array>
#include <cmath>

namespace utils
{
    std::array<double, 3> cross_product(const std::array<double, 3> &a, const std::array<double, 3> &b)
    {
        return {
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        };
    };

    double dot_product(const std::array<double, 3> &a, const std::array<double, 3> &b)
    {
        return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }

    std::array<double, 3> createVec(const std::array<double, 3> &a, const std::array<double, 3> &b)
    {
        return {
            b[0] - a[0],
            b[1] - a[1],
            b[2] - a[2],
        };
    }

    double magnitude(const std::array<double, 3> &a)
    {
        return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    }
}

#endif
