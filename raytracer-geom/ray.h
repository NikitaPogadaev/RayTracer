#pragma once

#include "vector.h"

class Ray {
public:
    Ray(const Vector& origin, const Vector& direction) : origin_(origin), direction_(direction) {
        direction_.Normalize();
    }
    Ray(const Ray&) = default;
    Ray(Ray&&) = default;
    Ray& operator=(const Ray&) = default;
    Ray& operator=(Ray&&) = default;
    ~Ray() = default;

    const Vector& GetOrigin() const {
        return origin_;
    }
    const Vector& GetDirection() const {
        return direction_;
    }

private:
    Vector origin_;
    Vector direction_;
};
