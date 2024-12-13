#pragma once

#include <array>
#include <cstddef>
#include <cmath>
#include <algorithm>

constexpr double kEpsilon = std::numeric_limits<double>::epsilon();

bool Eq(double aa, double bb) {
    return std::abs(aa - bb) < kEpsilon;
}

class Vector {
public:
    Vector() {
        data_[0] = 0.;
        data_[1] = 0.;
        data_[2] = 0.;
    }
    Vector(double x, double y, double z) {
        data_[0] = x;
        data_[1] = y;
        data_[2] = z;
    }
    Vector(const Vector&) = default;
    Vector(Vector&&) = default;
    Vector& operator=(const Vector&) = default;
    Vector& operator=(Vector&&) = default;
    ~Vector() = default;

    double& operator[](size_t ind) {
        return data_[ind];
    }
    double operator[](size_t ind) const {
        return data_[ind];
    }

    Vector operator-() const {
        return Vector(-data_[0], -data_[1], -data_[2]);
    }

    double Max() const {
        return std::max(std::max(data_[0], data_[1]), data_[2]);
    }

    double Min() const {
        return std::min(std::min(data_[0], data_[1]), data_[2]);
    }

    void Normalize() {
        double norm = std::sqrt(data_[0] * data_[0] + data_[1] * data_[1] + data_[2] * data_[2]);
        if (Eq(norm, 0.)) {
            return;
        }
        data_[0] /= norm;
        data_[1] /= norm;
        data_[2] /= norm;
    }

    bool operator==(const Vector& vec) const {
        return (Eq(data_[0], vec[0]) && Eq(data_[1], vec[1]) && Eq(data_[2], vec[2]));
    }

    Vector operator*(const Vector& vec) const {
        Vector tmp = vec;
        tmp[0] *= data_[0];
        tmp[1] *= data_[1];
        tmp[2] *= data_[2];
        return tmp;
    }

    Vector& operator*=(const Vector& vec) {
        data_[0] *= vec[0];
        data_[1] *= vec[1];
        data_[2] *= vec[2];
        return *this;
    }

    Vector operator/(const Vector& vec) const {
        Vector tmp = *this;
        tmp[0] /= vec[0];
        tmp[1] /= vec[1];
        tmp[2] /= vec[2];
        return tmp;
    }

    Vector& operator/=(const Vector& vec) {
        data_[0] /= vec[0];
        data_[1] /= vec[1];
        data_[2] /= vec[2];
        return *this;
    }

    Vector operator+(const Vector& vec) const {
        Vector tmp = vec;
        tmp[0] += data_[0];
        tmp[1] += data_[1];
        tmp[2] += data_[2];
        return tmp;
    }

    Vector& operator+=(const Vector& vec) {
        data_[0] += vec[0];
        data_[1] += vec[1];
        data_[2] += vec[2];
        return *this;
    }

    Vector operator-(const Vector& vec) const {
        Vector tmp = *this;
        tmp[0] -= vec[0];
        tmp[1] -= vec[1];
        tmp[2] -= vec[2];
        return tmp;
    }

    Vector& operator-=(const Vector& vec) {
        data_[0] -= vec[0];
        data_[1] -= vec[1];
        data_[2] -= vec[2];
        return *this;
    }

    Vector& operator*=(double lam) {
        data_[0] *= lam;
        data_[1] *= lam;
        data_[2] *= lam;
        return *this;
    }

    Vector& operator/=(double lam) {
        if (lam == 0.) {
            return *this;
        }
        data_[0] /= lam;
        data_[1] /= lam;
        data_[2] /= lam;
        return *this;
    }
    Vector operator*(double lam) const {
        Vector tmp = *this;
        tmp[0] *= lam;
        tmp[1] *= lam;
        tmp[2] *= lam;
        return tmp;
    }

    Vector operator/(double lam) const {
        Vector tmp = *this;
        if (lam == 0.) {
            return tmp;
        }
        tmp[0] /= lam;
        tmp[1] /= lam;
        tmp[2] /= lam;
        return tmp;
    }

private:
    std::array<double, 3> data_;
};

double DotProduct(const Vector& a, const Vector& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector CrossProduct(const Vector& a, const Vector& b) {
    return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}
double Length(const Vector& v) {
    return std::sqrt(DotProduct(v, v));
}
