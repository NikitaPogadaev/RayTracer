#pragma once

#include "material.h"
#include "vector.h"
#include "object.h"
#include "light.h"

#include <vector>
#include <unordered_map>
#include <string>
#include <sstream>
#include <filesystem>
#include <fstream>
#include <cctype>
#include <optional>
#include <iostream>
#include <deque>

struct VectorInds {
    int vertex;
    std::optional<int> texture;
    std::optional<int> normal;
};

std::vector<VectorInds> ParseFace(const std::string& fline) {
    std::vector<VectorInds> faceinds;
    std::istringstream ss(fline);
    std::string prefix;
    ss >> prefix;

    std::string element;
    while (ss >> element) {
        VectorInds fv;
        std::istringstream elstream(element);
        std::string vertex, texture, normal;
        std::getline(elstream, vertex, '/');
        int ind = std::stoi(vertex);
        if (ind > 0) {
            --ind;
        }
        fv.vertex = ind;

        if (std::getline(elstream, texture, '/')) {
            if (!texture.empty()) {
                ind = std::stoi(texture);
                if (ind > 0) {
                    --ind;
                }
                fv.texture = ind;
            }
            if (std::getline(elstream, normal, '/')) {
                if (!normal.empty()) {
                    ind = std::stoi(normal);
                    if (ind > 0) {
                        --ind;
                    }
                    fv.normal = ind;
                }
            }
        }
        faceinds.push_back(fv);
    }
    return faceinds;
}

class Scene {
public:
    const std::vector<Object>& GetObjects() const {
        return faces_;
    }
    const std::vector<SphereObject>& GetSphereObjects() const {
        return spheres_;
    }
    const std::vector<Light>& GetLights() const {
        return lights_;
    }
    const std::unordered_map<std::string, Material>& GetMaterials() const {
        return materials_;
    }

    friend Scene ReadScene(const std::filesystem::path& path);

private:
    std::deque<Vector> normals_;
    std::vector<Vector> verts_;
    std::vector<Object> faces_;
    std::vector<SphereObject> spheres_;
    std::vector<Light> lights_;
    std::unordered_map<std::string, Material> materials_;
};

std::unordered_map<std::string, Material> ReadMaterials(const std::filesystem::path& path) {
    std::ifstream file(path);
    std::unordered_map<std::string, Material> mtls;
    if (!file.is_open()) {
        return mtls;
    }
    std::string line;
    std::string name;
    Material mtr;
    mtr.albedo = Vector(1., 0., 0.);
    mtr.refraction_index = 1.;
    bool valid = false;
    while (std::getline(file, line)) {
        std::istringstream words(line);
        std::string token;
        while (words >> token) {
            if (token == "newmtl") {
                if (valid) {
                    mtr.name = name;
                    mtls[name] = mtr;
                } else {
                    valid = true;
                }
                mtr = Material();
                mtr.albedo = Vector(1., 0., 0.);
                mtr.refraction_index = 1.;
                words >> name;

                break;
            }
            if (token == "Ka" || token == "Ks" || token == "Kd" || token == "Ke" || token == "al") {
                std::string co1, co2, co3;
                words >> co1 >> co2 >> co3;
                Vector vec = Vector(std::stod(co1), std::stod(co2), std::stod(co3));
                if (token == "Ka") {
                    mtr.ambient_color = vec;
                } else if (token == "Ks") {
                    mtr.specular_color = vec;
                } else if (token == "Kd") {
                    mtr.diffuse_color = vec;
                } else if (token == "Ke") {
                    mtr.intensity = vec;
                } else {
                    mtr.albedo = vec;
                }
                break;
            }
            if (token == "Ns") {
                std::string co;
                words >> co;
                mtr.specular_exponent = std::stod(co);
                break;
            }
            if (token == "Ni") {
                std::string co;
                words >> co;
                mtr.refraction_index = std::stod(co);
                break;
            }
            break;
        }
    }
    if (valid) {
        mtr.name = name;
        mtls[name] = mtr;
    }
    file.close();
    return mtls;
}

Scene ReadScene(const std::filesystem::path& path) {
    std::ifstream file(path);
    Scene scene;
    if (!file.is_open()) {
        return scene;
    }
    std::string line;
    std::string name;
    while (std::getline(file, line)) {
        std::istringstream words(line);
        std::string token;
        while (words >> token) {
            if (token == "mtllib") {
                words >> token;
                std::filesystem::path mtlpath{token};
                scene.materials_ = ReadMaterials(path.parent_path() / mtlpath);
                break;
            }
            if (token == "v" || token == "vn") {
                std::string co1, co2, co3;
                words >> co1 >> co2 >> co3;
                if (token == "v") {
                    scene.verts_.emplace_back(std::stod(co1), std::stod(co2), std::stod(co3));
                } else {
                    scene.normals_.emplace_back(std::stod(co1), std::stod(co2), std::stod(co3));
                }
            }
            if (token == "usemtl") {
                words >> name;
                break;
            }
            if (token == "P") {
                std::string co1, co2, co3;
                std::string rgb1, rgb2, rgb3;
                words >> co1 >> co2 >> co3 >> rgb1 >> rgb2 >> rgb3;
                scene.lights_.emplace_back(
                    Vector(std::stod(co1), std::stod(co2), std::stod(co3)),
                    Vector(std::stod(rgb1), std::stod(rgb2), std::stod(rgb3)));
                break;
            }
            if (token == "S") {
                std::string co1, co2, co3;
                std::string rr;
                words >> co1 >> co2 >> co3 >> rr;
                scene.spheres_.emplace_back(
                    Sphere(Vector(std::stod(co1), std::stod(co2), std::stod(co3)), std::stod(rr)),
                    &(scene.materials_[name]));
                break;
            }
            if (token == "f") {
                std::vector<VectorInds> face = ParseFace(line);
                if (face.size() < 3) {
                    break;
                }
                VectorInds firstvert = face[0];
                VectorInds secondvert = face[1];
                int len = scene.verts_.size();
                int lenlen = scene.normals_.size();
                for (size_t i = 2; i < face.size(); ++i) {
                    Triangle tr(scene.verts_[((firstvert.vertex % len) + len) % len],
                                scene.verts_[((secondvert.vertex % len) + len) % len],
                                scene.verts_[((face[i].vertex % len) + len) % len]);
                    const Material* mtr = &(scene.materials_[name]);
                    Vector* ptr1 =
                        (firstvert.normal.has_value()
                             ? &(scene.normals_[(((firstvert.normal.value() % lenlen) + lenlen) %
                                                 lenlen)])
                             : nullptr);
                    Vector* ptr2 =
                        (secondvert.normal.has_value()
                             ? &(scene.normals_[(((secondvert.normal.value() % lenlen) + lenlen) %
                                                 lenlen)])
                             : nullptr);
                    Vector* ptr3 =
                        (face[i].normal.has_value()
                             ? &(scene.normals_[(((face[i].normal.value() % lenlen) + lenlen) %
                                                 lenlen)])
                             : nullptr);
                    scene.faces_.emplace_back(tr, mtr, ptr1, ptr2, ptr3);
                    secondvert = face[i];
                }
            }
            break;
        }
    }
    file.close();
    return scene;
}
