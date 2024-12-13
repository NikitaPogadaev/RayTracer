#pragma once

#include "vector.h"

struct Light {
    Light(const Vector& pos, const Vector& in) : position(pos), intensity(in) {
    }
    Vector position;
    Vector intensity;
    Light(const Light&) = default;
    Light(Light&&) = default;
    Light& operator=(const Light&) = default;
    Light& operator=(Light&&) = default;
    ~Light() = default;
};
