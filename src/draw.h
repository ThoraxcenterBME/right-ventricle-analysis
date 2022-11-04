#pragma once
// Disable warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/glm.hpp>
DISABLE_WARNINGS_POP()
#include <framework/mesh.h>
#include <span>

void drawMeshWithColors(const Mesh& mesh, std::span<const glm::vec3> colors);