#pragma once
#include "image.h"
// Suppress warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
DISABLE_WARNINGS_POP()
#include <filesystem>
#include <optional>
#include <span>
#include <vector>
#include <map>

struct Vertex {
	glm::vec3 position;
	glm::vec3 normal;
	glm::vec2 texCoord; // Texture coordinate
    std::vector<int> ring; // Saves the vertices that make the 1-ring neighbourhood
	int ringCount; 
	double curvature;
    double indexed_curv;
	int index; 
	void setCurvature(double c) {
		curvature = c; 
	}
	void set_index_curv(double c) {
		indexed_curv = c; 
	}
    bool exclude = false; // whether it should be excluded in the curvature calculation
	[[nodiscard]] constexpr bool operator==(const Vertex&) const noexcept = default;
};

struct Material {
	glm::vec3 kd; // Diffuse color.
	glm::vec3 ks{ 0.0f };
	float shininess{ 1.0f };
	float transparency{ 1.0f };

	// Optional texture that replaces kd; use as follows:
	// 
	// if (material.kdTexture) {
	//   material.kdTexture->getTexel(...);
	// }
	std::shared_ptr<Image> kdTexture;
};

struct Mesh {
	// Vertices contain the vertex positions and normals of the mesh.
	std::vector<Vertex> vertices;
	// A triangle contains a triplet of values corresponding to the indices of the 3 vertices in the vertices array.
	std::vector<glm::uvec3> triangles;
	Material material;
    std::map<int, std::vector<int>> vertexToTri;
    double radius; 
};

[[nodiscard]] std::vector<Mesh> loadMesh(const std::filesystem::path& file, bool normalize = false);
[[nodiscard]] Mesh mergeMeshes(std::span<const Mesh> meshes);
void meshFlipX(Mesh& mesh);
void meshFlipY(Mesh& mesh);
void meshFlipZ(Mesh& mesh);