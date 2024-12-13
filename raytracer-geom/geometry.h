#pragma once

#include "vector.h"
#include "sphere.h"
#include "intersection.h"
#include "triangle.h"
#include "ray.h"

#include <optional>

enum class RootCount { kZero, kOne, kTwo, kInf };

struct Roots {
    RootCount count;
    double first;
    double second;
};

Roots SolveQuadratic(double a, double b, double c) {
    if (Eq(a, 0.)) {
        if (Eq(b, 0.)) {
            return {(Eq(c, 0.) ? RootCount::kInf : RootCount::kZero), 0., 0.};
        }
        return {RootCount::kOne, -c / b, -c / b};
    }

    if (b * b < 4. * a * c) {
        return {RootCount::kZero, 0., 0.};
    }
    double d = b * b - 4. * a * c;
    double d1 = (-b - sqrt(d)) / (2. * a);
    double d2 = (-b + sqrt(d)) / (2. * a);
    if (d1 > d2) {
        std::swap(d1, d2);
    }

    return {RootCount::kTwo, d1, d2};
}

std::optional<Intersection> GetIntersection(const Ray& ray, const Sphere& sphere) {
    Vector o_pvec = ray.GetOrigin() - sphere.GetCenter();
    double op = Length(o_pvec);
    double aq = DotProduct(ray.GetDirection(), ray.GetDirection());
    double bq = 2. * DotProduct(o_pvec, ray.GetDirection());
    double cq = DotProduct(o_pvec, o_pvec) - sphere.GetRadius() * sphere.GetRadius();
    auto root = SolveQuadratic(aq, bq, cq);
    if (root.count == RootCount::kZero || root.count == RootCount::kInf) {
        return std::nullopt;
    }

    if (root.first >= 0.) {
        Vector peres = ray.GetOrigin() + ray.GetDirection() * root.first;
        Vector norm = peres - sphere.GetCenter();
        double rast = Length(peres - ray.GetOrigin());
        if (op <= sphere.GetRadius()) {
            norm = -norm;
        }
        return Intersection(peres, norm, rast);
    }
    if (root.second >= 0.) {
        Vector peres = ray.GetOrigin() + ray.GetDirection() * root.second;
        Vector norm = peres - sphere.GetCenter();
        double rast = Length(peres - ray.GetOrigin());
        if (op <= sphere.GetRadius()) {
            norm = -norm;
        }
        return Intersection(peres, norm, rast);
    }
    return std::nullopt;
}

std::optional<Intersection> GetIntersection(const Ray& ray, const Triangle& triangle) {
    auto ray_vector = ray.GetDirection();
    auto ray_origin = ray.GetOrigin();
    Vector edge1 = triangle[1] - triangle[0];
    Vector edge2 = triangle[2] - triangle[0];
    Vector ray_cross_e2 = CrossProduct(ray_vector, edge2);
    double det = DotProduct(edge1, ray_cross_e2);

    if (det > -kEpsilon && det < kEpsilon) {
        return std::nullopt;
    }

    double inv_det = 1.0 / det;
    Vector s = ray_origin - triangle[0];
    double u = inv_det * DotProduct(s, ray_cross_e2);

    if (u < 0. || u > 1.) {
        return std::nullopt;
    }

    Vector s_cross_e1 = CrossProduct(s, edge1);
    double v = inv_det * DotProduct(ray_vector, s_cross_e1);

    if (v < 0. || u + v > 1) {
        return std::nullopt;
    }

    double t = inv_det * DotProduct(edge2, s_cross_e1);

    if (t > kEpsilon) {
        auto crans = CrossProduct(edge1, edge2);
        auto ans = ray_origin + ray_vector * t;
        if (DotProduct(ray_vector, crans) > 0.) {
            crans = -crans;
        }
        return Intersection(ans, crans, Length(ray_origin - ans));
    } else {
        return std::nullopt;
    }
}

Vector Reflect(const Vector& ray, const Vector& normal) {
    Vector diff = normal * (2. * DotProduct(ray, normal));
    return ray - diff;
}
std::optional<Vector> Refract(const Vector& ray, const Vector& normal, double eta) {
    if (Length(CrossProduct(ray, normal)) * eta >= 1.) {
        return std::nullopt;
    }
    double cso = -DotProduct(normal, ray);
    return ray * eta + normal * (eta * cso - std::sqrt(1. - eta * eta * (1. - cso * cso)));
}
Vector GetBarycentricCoords(const Triangle& triangle, const Vector& point) {
    double hole = triangle.Area();
    double barx = Triangle(triangle[1], triangle[2], point).Area() / hole;
    double bary = Triangle(triangle[0], triangle[2], point).Area() / hole;
    double barz = Triangle(triangle[1], triangle[0], point).Area() / hole;
    return Vector(barx, bary, barz);
}
