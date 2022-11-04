#pragma once
// Disable warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/geometric.hpp>
#include <glm/glm.hpp>
DISABLE_WARNINGS_POP()
#include <iostream>
#include <limits>
#include <utility>

using Color = glm::vec3;

// !!! DO NOT MODIFY THIS STRUCT !!!
struct MaterialInformation {
    Color Kd { 0.5f, 0.5f, 0.5f }; // Diffuse coefficient per vertex.
    Color Ks { 0.5f, 0.5f, 0.5f }; // Specularity coefficient per vertex.
    float shininess = 20.0f; // Exponent for Phong and Blinn-Phong specularities per vertex.

    // Gooch shading
    float goochB, goochY, goochAlpha, goochBeta;
};
struct Positions {
    glm::vec3 vertex; // Vertex position in world space
    glm::vec3 light; // Light position in world space
    glm::vec3 camera; // Camera position in world space
};

Color debugColor(const MaterialInformation& materialInformation, const Positions& positions, const glm::vec3& normal, const Color& lightColor)
{
    // This function you can use in any way you like!
    // E.g., for debugging purposes!
    return (normal + 1.0f) / 2.0f;

    // or random color per vertex:
    // const uint32_t hashX = std::hash<float>()(positions.vertex.x);
    // const uint32_t hashY = std::hash<float>()(positions.vertex.y);
    // const uint32_t hashZ = std::hash<float>()(positions.vertex.z);
    // return Color {
    //     (double)hashX / std::numeric_limits<uint32_t>::max(),
    //     (double)hashY / std::numeric_limits<uint32_t>::max(),
    //     (double)hashZ / std::numeric_limits<uint32_t>::max()
    // };

    // or material information:
    // return materialInformation.Kd;
}

// Standard lambertian shading: Kd * dot(N,L), clamped to zero when negative. Where L is the light vector.
Color diffuseOnly(const MaterialInformation& materialInformation, const Positions& positions, const glm::vec3& normal, const Color& lightColor)
{
    return glm::vec3(0, 0.5, 0.5);
}

Color phongSpecularOnly(const MaterialInformation& materialInformation, const Positions& positions, const glm::vec3& normal, const Color& lightColor)
{
    return glm::vec3(0, 1, 0);
}

Color blinnPhongSpecularOnly(const MaterialInformation& materialInformation, const Positions& positions, const glm::vec3& normal, const Color& lightColor)
{
    return glm::vec3(0, 0, 1);
}

Color gooch(const MaterialInformation& materialInformation, const Positions& positions, const glm::vec3& normal, const Color& lightColor)
{
    return glm::vec3(0, 0, 1);
}

glm::vec3 userInteractionSphere(const glm::vec3& selectedPosition, const glm::vec3& cameraPosition)
{
    return glm::vec3(1, 1, 1);
}
