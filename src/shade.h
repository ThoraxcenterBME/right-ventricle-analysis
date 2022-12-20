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

// Utility function for calculating a reflection vector
// The incident vector points _towards_ the light point
// The normal vector _must_ be normalized
glm::vec3 calculateReflection(const glm::vec3& incident, const glm::vec3& normal)
{
    glm::vec3 reflect = -1.0f * incident + 2.0f * glm::dot(incident, normal) * normal;

    return reflect;
}


// Standard lambertian shading: Kd * dot(N,L), clamped to zero when negative. Where L is the light vector.
glm::vec3 diffuseOnly(const MaterialInformation& shadingData, const Positions& positions, const glm::vec3& normal, const glm::vec3& lightColor)
{
    // Calculate the diffuse coefficient per vertex
    glm::vec3 light_dir = glm::normalize(positions.light - positions.vertex);

    float x_diffuse = shadingData.Kd[0] * glm::max(0.0f, glm::dot(light_dir, glm::normalize(normal)));
    float y_diffuse = shadingData.Kd[1] * glm::max(0.0f, glm::dot(light_dir, glm::normalize(normal)));
    float z_diffuse = shadingData.Kd[2] * glm::max(0.0f, glm::dot(light_dir, glm::normalize(normal)));

    return glm::vec3(x_diffuse, y_diffuse, z_diffuse);
}

// Phong (!) Shading Specularity (http://en.wikipedia.org/wiki/Blinn%E2%80%93Phong_shading_model)
// Follow the course, only calculate Ks pow(dot(V,R),shininess), where V is the view vector and R is the Reflection vector of the light (like in pool billard computed from the LightPos, vertexPos and normal).
// When computing specularities like this, verify that the light is on the right side of the surface, with respect to the normal
// E.g., for a plane, the light source below the plane cannot cast light on the top, hence, there can also not be any specularity.
glm::vec3 phongSpecularOnly(const MaterialInformation& shadingData, const glm::vec3& vertexPos, const glm::vec3& normal, const glm::vec3& lightPos, const glm::vec3& cameraPos)
{
    // Check if light source is above the vertex. Use cosine and normalized.
    glm::vec3 light = glm::normalize(glm::vec3(lightPos[0], lightPos[1], lightPos[2])
        - glm::vec3(vertexPos[0], vertexPos[1], vertexPos[2]));

    if (glm::dot(light, glm::normalize(normal)) <= 0) {
        return glm::vec3(0, 0, 0);
    }

    // Find the reflection vector
    glm::vec3 reflect = glm::normalize(calculateReflection(light, glm::normalize(normal)));
    // Find the viewing vector
    glm::vec3 view = glm::normalize(glm::vec3(cameraPos[0], cameraPos[1], cameraPos[2])
        - glm::vec3(vertexPos[0], vertexPos[1], vertexPos[2]));

    // Calculate for each point, using the colours
    float x_spec = shadingData.Ks[0] * glm::pow(glm::max(0.0f, glm::dot(view, reflect)), shadingData.shininess);
    float y_spec = shadingData.Ks[1] * glm::pow(glm::max(0.0f, glm::dot(view, reflect)), shadingData.shininess);
    float z_spec = shadingData.Ks[2] * glm::pow(glm::max(0.0f, glm::dot(view, reflect)), shadingData.shininess);

    return glm::vec3(x_spec, y_spec, z_spec);
}

// Blinn-Phong Shading Specularity (http://en.wikipedia.org/wiki/Blinn%E2%80%93Phong_shading_model)
// Be careful!!! The pseudo code does some additional modifications to the formula seen in the course.
// Follow the course version and calculate ONLY Ks * pow(dot(N,H), shininess). The definition of H is given on the page above and in the course.
// The same test as before should be used.
glm::vec3 blinnPhongSpecularOnly(const MaterialInformation& shadingData, const glm::vec3& vertexPos, const glm::vec3& normal, const glm::vec3& lightPos, const glm::vec3& cameraPos)
{
    // Check if light source is above the vertex. Use cosine and normalized.
    glm::vec3 light = glm::normalize(glm::vec3(lightPos[0], lightPos[1], lightPos[2])
        - glm::vec3(vertexPos[0], vertexPos[1], vertexPos[2]));

    if (glm::dot(light, glm::normalize(normal)) <= 0) {
        return glm::vec3(0, 0, 0);
    }

    // Find the view vector
    glm::vec3 view = glm::normalize(glm::vec3(cameraPos[0], cameraPos[1], cameraPos[2])
        - glm::vec3(vertexPos[0], vertexPos[1], vertexPos[2]));
    // Find H vector - Blinn-Phong
    glm::vec3 h = glm::normalize(view + light);

    float x_bspec = shadingData.Ks[0] * glm::pow(glm::max(0.0f, glm::dot(glm::normalize(normal), h)), shadingData.shininess);
    float y_bspec = shadingData.Ks[1] * glm::pow(glm::max(0.0f, glm::dot(glm::normalize(normal), h)), shadingData.shininess);
    float z_bspec = shadingData.Ks[2] * glm::pow(glm::max(0.0f, glm::dot(glm::normalize(normal), h)), shadingData.shininess);

    return glm::vec3(x_bspec, y_bspec, z_bspec);
}

float obtainColour(const float kd_value, const int n, const float diffuse)
{
    // Special case for extreme value
    if (diffuse == 0.0f)
        return 0.0f;

    // Stores the possible colour values
    float d = static_cast<float>(kd_value / n);

    for (int i = 0; i < n; i++) {
        if (diffuse >= d * i && diffuse < d * (i + 1)) {
            return (d * i + d * (i + 1)) / 2;
        }
    }

    return kd_value; // otherwise, just return Kd
}


Color gooch(const MaterialInformation& materialInformation, const Positions& positions, const glm::vec3& normal, const Color& lightColor)
{
    return glm::vec3(0, 0, 1);
}


// This function finds the values t_1 and t_2 that will allow the ray to intersect the sphere.
std::vector<float> findOffsets(const glm::vec3& direction, const glm::vec3& init)
{
    // Rewrite in form a * t^2 + b * t + (c - 1.5^2) = 0
    float a = glm::dot(direction, direction);
    float b = 2.0f * glm::dot(direction, init);
    float c = glm::dot(init, init) - (1.5f * 1.5f);

    // Define sqrt(b^2 - 4ac)
    float sqrt_discr = glm::sqrt((b * b) - (4.0f * a * c));

    // Define (-b +- discr) / (2a)
    float t_1 = static_cast<float>((-b + sqrt_discr) / (2.0f * a));
    float t_2 = static_cast<float>((-b - sqrt_discr) / (2.0f * a));

    return std::vector<float> { t_1, t_2 };
}

// This function determine which of the intersection points to use. Closer to the camera position.
glm::vec3 findCloser(const glm::vec3& direction, const glm::vec3& init, const std::vector<float>& offsets)
{
    float t_1 = offsets[0];
    float t_2 = offsets[1];

    glm::vec3 point_1 = init + glm::vec3(t_1 * direction.x, t_1 * direction.y, t_1 * direction.z);
    glm::vec3 point_2 = init + glm::vec3(t_2 * direction.x, t_2 * direction.y, t_2 * direction.z);

    // Find which one is closer to init in absolute distance
    float distance_1 = glm::distance(init, point_1);
    float distance_2 = glm::distance(init, point_2);

    // Return the closer point
    if (distance_1 < distance_2) {
        return point_1;
    } else {
        return point_2;
    }
}

// RETURN the new light position, defined as follows:
// selectedPos is a location on the mesh. Use this location to place the light source to cover the location as seen from camPos (cover the cursor).
// Further, the light should be at a distance of 1.5 from the origin of the scene - in other words, located on a sphere of radius 1.5 around the origin.
// The selectedPos is guaranteed to always be within the sphere.
glm::vec3 userInteractionSphere(const glm::vec3& selectedPos, const glm::vec3& camPos)
{
    glm::vec3 direction = selectedPos - camPos;
    std::vector<float> offsets = findOffsets(direction, camPos);
    glm::vec3 light_pos = findCloser(direction, camPos, offsets);

    return light_pos;
}