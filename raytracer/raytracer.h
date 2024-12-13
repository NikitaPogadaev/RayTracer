#pragma once

#include "options/camera_options.h"
#include "options/render_options.h"
#include "image.h"
#include "scene.h"
#include "geometry.h"

#include <filesystem>

bool IsVisible(const Light& light, Scene& scene, Vector position) {
    if (position == light.position) {
        return true;
    }
    Ray dir(light.position, position - light.position);
    double dist = Length(position - light.position);
    for (const auto& obj : scene.GetObjects()) {
        auto inter = GetIntersection(dir, obj.polygon);
        if (inter == std::nullopt) {
            continue;
        }
        if (inter.value().GetDistance() < dist) {
            return false;
        }
    }
    for (const auto& obj : scene.GetSphereObjects()) {
        auto inter = GetIntersection(dir, obj.sphere);
        if (inter == std::nullopt) {
            continue;
        }
        if (inter.value().GetDistance() < dist) {
            return false;
        }
    }
    return true;
}

Vector ComputePreRGB(Scene& scene, Ray& ray, int depth) {
    Vector pixel(0., 0., 0.);
    if (scene.GetObjects().empty() && scene.GetSphereObjects().empty()) {
        return pixel;
    }
    const Object* pptr = nullptr;
    const SphereObject* sptr = nullptr;
    double pdist = -1.;
    double sdist = -1.;
    Vector ppos;
    Vector spos;
    Vector pnorm;
    Vector snorm;
    for (const auto& obj : scene.GetObjects()) {
        auto inter = GetIntersection(ray, obj.polygon);
        if (inter == std::nullopt) {
            continue;
        } else {
            double dist = inter.value().GetDistance();
            if (pptr == nullptr || dist < pdist) {
                pdist = dist;
                pptr = &obj;
                ppos = inter.value().GetPosition();
                pnorm = inter.value().GetNormal();
            }
        }
    }
    for (const auto& obj : scene.GetSphereObjects()) {
        auto inter = GetIntersection(ray, obj.sphere);
        if (inter == std::nullopt) {
            continue;
        } else {
            double dist = inter.value().GetDistance();
            if (sptr == nullptr || dist < sdist) {
                sdist = dist;
                sptr = &obj;
                spos = inter.value().GetPosition();
                snorm = inter.value().GetNormal();
            }
        }
    }

    if (sdist == -1. && pdist == -1.) {
        return pixel;
    }

    if (pptr != nullptr && pptr->GetNormal(0) != nullptr) {
        Vector coords = GetBarycentricCoords(pptr->polygon, ppos);
        pnorm = *(pptr->GetNormal(0)) * coords[0] + *(pptr->GetNormal(1)) * coords[1] +
                *(pptr->GetNormal(2)) * coords[2];
    }

    if (pptr != nullptr && (sdist == -1. || pdist < sdist)) {
        pixel = pixel + pptr->material->ambient_color + pptr->material->intensity;
        Vector shiftppos = ppos + pnorm * 0.00001;
        Vector antippos = ppos - pnorm * 0.00001;
        for (const auto& light : scene.GetLights()) {
            if (!IsVisible(light, scene, shiftppos)) {
                continue;
            }
            Vector lightdir = light.position - shiftppos;
            lightdir.Normalize();
            Vector vr = Reflect(-lightdir, pnorm);
            vr.Normalize();
            if (pptr->material->albedo[0] > 0.) {
                pixel = pixel + (light.intensity * pptr->material->diffuse_color *
                                     std::max(0., DotProduct(pnorm, lightdir)) +
                                 light.intensity *
                                     std::pow(std::max(0., DotProduct(-ray.GetDirection(), vr)),
                                              pptr->material->specular_exponent) *
                                     pptr->material->specular_color) *
                                    pptr->material->albedo[0];
            }
        }                 ////////////////////
        if (depth > 1) {  //???
            if (pptr->material->albedo[1] > 0.) {
                Ray vrray(shiftppos, Reflect(ray.GetDirection(), pnorm));
                pixel = pixel + ComputePreRGB(scene, vrray, depth - 1) * pptr->material->albedo[1];
            }
            if (pptr->material->albedo[2] > 0.) {
                auto ref =
                    Refract(ray.GetDirection(), pnorm, 1. / pptr->material->refraction_index);
                if (ref != std::nullopt) {
                    Ray refray(antippos, ref.value());
                    pixel =
                        pixel + ComputePreRGB(scene, refray, depth - 1) * pptr->material->albedo[2];
                }
            }
        }
    } else {
        pixel = pixel + sptr->material->ambient_color + sptr->material->intensity;
        Vector shiftspos = spos + snorm * 0.00001;
        Vector antispos = spos - snorm * 0.00001;
        for (const auto& light : scene.GetLights()) {
            if (!IsVisible(light, scene, shiftspos)) {
                continue;
            }
            Vector lightdir = light.position - shiftspos;
            lightdir.Normalize();
            Vector vr = Reflect(-lightdir, snorm);
            vr.Normalize();
            if (sptr->material->albedo[0] > 0.) {
                pixel = pixel + (light.intensity * sptr->material->diffuse_color *
                                     std::max(0., DotProduct(snorm, lightdir)) +
                                 light.intensity *
                                     std::pow(std::max(0., DotProduct(-ray.GetDirection(), vr)),
                                              sptr->material->specular_exponent) *
                                     sptr->material->specular_color) *
                                    sptr->material->albedo[0];
            }
        }
        if (depth > 1) {  //???
            bool isin = (DotProduct(spos - sptr->sphere.GetCenter(), spos - snorm) > 0.);
            if (isin) {
                if (sptr->material->albedo[2] > 0.) {
                    auto ref = Refract(ray.GetDirection(), snorm, sptr->material->refraction_index);
                    if (ref != std::nullopt) {
                        Ray refray(antispos, ref.value());
                        pixel = pixel + ComputePreRGB(scene, refray, depth - 1);
                    }
                }
            } else {
                if (sptr->material->albedo[1] > 0.) {
                    Ray vrray(shiftspos, Reflect(ray.GetDirection(), snorm));
                    pixel =
                        pixel + ComputePreRGB(scene, vrray, depth - 1) * sptr->material->albedo[1];
                }
                if (sptr->material->albedo[2] > 0.) {
                    auto ref =
                        Refract(ray.GetDirection(), snorm, 1. / sptr->material->refraction_index);
                    if (ref != std::nullopt) {
                        Ray refray(antispos, ref.value());
                        pixel = pixel +
                                ComputePreRGB(scene, refray, depth - 1) * sptr->material->albedo[2];
                    }
                }
            }
        }
    }
    return pixel;
}

Image Render(const std::filesystem::path& path, const CameraOptions& camera_options,
             const RenderOptions& render_options) {
    Image img(camera_options.screen_width, camera_options.screen_height);
    Scene scene = ReadScene(path);
    std::vector<std::vector<Vector>> screen(camera_options.screen_height,
                                            std::vector<Vector>(camera_options.screen_width));

    double height = camera_options.screen_height;
    double width = camera_options.screen_width;
    // double screendist = height / 2. / std::tan(camera_options.fov / 2.);
    double hcoef = 2. * std::tan(camera_options.fov / 2.);
    double wcoef = hcoef * width / height;

    Vector forward = camera_options.look_to - camera_options.look_from;
    forward.Normalize();
    Vector right(1., 0, 0);
    Vector up(0., 1., 0.);
    if (Eq(Length(CrossProduct(up, forward)), 0.)) {
        if (Eq(up[1], forward[1])) {
            right = Vector(1., 0., 0.);
            up = Vector(0., 0., -1.);
            forward = Vector(0., 1., 0.);
        } else {
            right = Vector(1., 0., 0.);
            up = Vector(0., 0., 1.);
            forward = Vector(0., -1., 0.);
        }
    } else {
        right = CrossProduct(up, forward);
        right.Normalize();
        up = CrossProduct(right, forward);
    }

    if (render_options.mode == RenderMode::kDepth) {
        double maxdist = -1.;
        std::vector<std::vector<bool>> screenflag(
            camera_options.screen_height, std::vector<bool>(camera_options.screen_width, false));
        for (int j = 0; j < camera_options.screen_height; ++j) {
            for (int i = 0; i < camera_options.screen_width; ++i) {
                screen[j][i] = {1., 1., 1.};
                double x = (i + 0.5 - width / 2.0) * wcoef / width;
                double y = (j + 0.5 - height / 2.0) * hcoef / height;

                Vector raydir = (forward + right * x + up * y);
                Ray nextray(camera_options.look_from, raydir);
                double dist, mindist = -1.;
                for (const auto& obj : scene.GetObjects()) {
                    auto inter = GetIntersection(nextray, obj.polygon);
                    if (inter == std::nullopt) {
                        continue;
                    } else {
                        dist = inter.value().GetDistance();
                        if (mindist == -1. || dist < mindist) {
                            mindist = dist;
                        }
                    }
                }

                for (const auto& obj : scene.GetSphereObjects()) {
                    auto inter = GetIntersection(nextray, obj.sphere);
                    if (inter == std::nullopt) {
                        continue;
                    } else {
                        dist = inter.value().GetDistance();
                        if (mindist == -1. || dist < mindist) {
                            mindist = dist;
                        }
                    }
                }
                maxdist = std::max(mindist, maxdist);
                if (mindist > 0) {
                    screen[j][i] = {mindist, mindist, mindist};
                    screenflag[j][i] = true;
                }
            }
        }

        for (int j = 0; j < camera_options.screen_height; ++j) {
            for (int i = 0; i < camera_options.screen_width; ++i) {
                double ncoef = screenflag[j][i] ? maxdist : 1.;
                RGB pixel(std::lround(screen[j][i][0] / ncoef * 255.),
                          std::lround(screen[j][i][0] / ncoef * 255.),
                          std::lround(screen[j][i][0] / ncoef * 255.));
                img.SetPixel(pixel, j, camera_options.screen_width - 1 - i);
            }
        }
        return img;

    } else if (render_options.mode == RenderMode::kNormal) {
        for (int j = 0; j < camera_options.screen_height; ++j) {
            for (int i = 0; i < camera_options.screen_width; ++i) {
                screen[j][i] = {0., 0., 0.};
                double x = (i + 0.5 - width / 2.0) * wcoef / width;
                double y = (j + 0.5 - height / 2.0) * hcoef / height;

                Vector raydir = (forward + right * x + up * y);
                Ray nextray(camera_options.look_from, raydir);
                double dist, mindist = -1.;
                for (const auto& obj : scene.GetObjects()) {
                    auto inter = GetIntersection(nextray, obj.polygon);
                    if (inter == std::nullopt) {
                        continue;
                    } else {
                        dist = inter.value().GetDistance();
                        if (mindist == -1. || dist < mindist) {
                            mindist = dist;
                            if (obj.GetNormal(0) == nullptr) {
                                screen[j][i] =
                                    inter.value().GetNormal() * 0.5 + Vector(0.5, 0.5, 0.5);
                            } else {
                                Vector coords =
                                    GetBarycentricCoords(obj.polygon, inter.value().GetPosition());
                                coords = *obj.GetNormal(0) * coords[0] +
                                         *obj.GetNormal(1) * coords[1] +
                                         *obj.GetNormal(2) * coords[2];
                                screen[j][i] = coords * 0.5 + Vector(0.5, 0.5, 0.5);
                            }
                        }
                    }
                }

                for (const auto& obj : scene.GetSphereObjects()) {
                    auto inter = GetIntersection(nextray, obj.sphere);
                    if (inter == std::nullopt) {
                        continue;
                    } else {
                        dist = inter.value().GetDistance();
                        if (mindist == -1. || dist < mindist) {
                            mindist = dist;
                            screen[j][i] = inter.value().GetNormal() * 0.5 + Vector(0.5, 0.5, 0.5);
                        }
                    }
                }
                RGB pixel(std::lround(screen[j][i][0] * 255.), std::lround(screen[j][i][1] * 255.),
                          std::lround(screen[j][i][2] * 255.));
                img.SetPixel(pixel, j, camera_options.screen_width - 1 - i);
            }
        }
        return img;
    }

    ////////////////////////////////////////
    double ncoef = std::numeric_limits<double>::min();
    for (int j = 0; j < camera_options.screen_height; ++j) {
        for (int i = 0; i < camera_options.screen_width; ++i) {
            screen[j][i] = {0., 0., 0.};
            double x = (i + 0.5 - width / 2.0) * wcoef / width;
            double y = (j + 0.5 - height / 2.0) * hcoef / height;

            Vector raydir = (forward + right * x + up * y);
            Ray nextray(camera_options.look_from, raydir);
            screen[j][i] = ComputePreRGB(scene, nextray, render_options.depth);
            ncoef = std::max(ncoef, screen[j][i].Max());
        }
    }

    for (int j = 0; j < camera_options.screen_height; ++j) {
        for (int i = 0; i < camera_options.screen_width; ++i) {
            Vector tmp = screen[j][i];
            tmp = tmp * (tmp / (ncoef * ncoef) + Vector(1., 1., 1.)) / (tmp + Vector(1., 1., 1.));
            double ex = 1. / 2.2;
            tmp = Vector(std::pow(tmp[0], ex), std::pow(tmp[1], ex), std::pow(tmp[2], ex));
            RGB pixel(std::lround(tmp[0] * 255.), std::lround(tmp[1] * 255.),
                      std::lround(tmp[2] * 255.));
            img.SetPixel(pixel, j, camera_options.screen_width - 1 - i);
        }
    }
    return img;
}
