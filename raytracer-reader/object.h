#pragma once

#include "triangle.h"
#include "material.h"
#include "sphere.h"
#include "vector.h"

struct Object {
    Object(const Triangle& tr, const Material* mtr, const Vector* norm1 = nullptr,
           const Vector* norm2 = nullptr, const Vector* norm3 = nullptr)
        : material(mtr), polygon(tr) {
        normals_[0] = norm1;
        normals_[1] = norm2;
        normals_[2] = norm3;
    }

    // Object(const Object&) = default;
    // Object(Object&&) = default;
    // Object& operator=(const Object&) = default;
    // Object& operator=(Object&&) = default;
    // ~Object() = default;

    const Material* material = nullptr;
    Triangle polygon;

    const Vector* GetNormal(size_t index) const {
        if (index > 2) {
            return nullptr;
        }
        return normals_[index];
    }

private:
    std::array<const Vector*, 3> normals_;
};

struct SphereObject {
    SphereObject(const Sphere& sphere, const Material* mtl) : material(mtl), sphere(sphere) {
    }

    // SphereObject(const SphereObject&) = default;
    // SphereObject(SphereObject&&) = default;
    // SphereObject& operator=(const SphereObject&) = default;
    // SphereObject& operator=(SphereObject&&) = default;
    // ~SphereObject() = default;

    const Material* material = nullptr;
    Sphere sphere;
};
