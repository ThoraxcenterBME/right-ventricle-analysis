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

// Enumerated type for indicating the region the vertex belongs to 
enum Region {
	// IFW = Interior Free Wall
	// LFW = Laterior Free Wall
	// AFW = Anterior Free Wall 
	// SP = Septal body
	// IGN = Ignored 
	// UD = Initial value, so we haven't classified it 
	IFW=0, LFW=1, AFW=2, SP=3, IGN=4, UD=5
};

struct Vertex {
	// Variables for vertex 
	glm::vec3 position;
	glm::vec3 normal;
	glm::vec2 texCoord; // Texture coordinate
    std::vector<int> ring; // Saves the vertices that make the 1-ring neighbourhood
	int ringCount; 
	double curvature;
    double indexed_curv;
    bool exclude = false; // whether it should be excluded in the curvature calculation
    [[nodiscard]] constexpr bool operator==(const Vertex&) const noexcept = default;
	int index; 
    Region region = Region::UD;  
	
	// Setters (might remove)
	void setCurvature(double c) {
		curvature = c; 
	}
	void set_index_curv(double c) {
		indexed_curv = c; 
	}
	void set_region(Region r) {
		region = r; 
	}
	double get_indexed_curvature() {
		return indexed_curv; 
	}
};

struct Material {
	glm::vec3 kd; // Diffuse color.
	glm::vec3 ks{ 0.0f };
	float shininess{ 1.0f };
	float transparency{ 1.0f };
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
