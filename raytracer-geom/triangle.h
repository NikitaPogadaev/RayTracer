#pragma once

#include "vector.h"

#include <cstddef>
#include <algorithm>

class Triangle {
public:
    Triangle(const Vector& a, const Vector& b, const Vector& c) {
        data_[0] = a;
        data_[1] = b;
        data_[2] = c;
    }

    Triangle(const Triangle&) = default;
    Triangle(Triangle&&) = default;
    Triangle& operator=(const Triangle&) = default;
    Triangle& operator=(Triangle&&) = default;
    ~Triangle() = default;

    const Vector& operator[](size_t ind) const {
        return data_[ind];
    }
    double Area() const {
        return Length(CrossProduct(data_[1] - data_[0], data_[2] - data_[0])) / 2.;
    }

private:
    std::array<Vector, 3> data_;
};
