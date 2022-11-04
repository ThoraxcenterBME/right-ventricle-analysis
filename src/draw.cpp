#include "draw.h"
// Disable warnings in third-party code.
#include <framework/disable_all_warnings.h>
DISABLE_WARNINGS_PUSH()
#include <glm/gtc/type_ptr.hpp>
DISABLE_WARNINGS_POP()
#include <framework/window.h>

void drawMeshWithColors(const Mesh& mesh, std::span<const glm::vec3> colors)
{
    glBegin(GL_TRIANGLES);
    for (const auto& triangle : mesh.triangles) {
        for (int v = 0; v < 3; v++) {
            const auto& color = colors[triangle[v]];
            const auto& vertex = mesh.vertices[triangle[v]];

            glColor3fv(glm::value_ptr(color));
            glNormal3fv(glm::value_ptr(vertex.normal));
            glVertex3fv(glm::value_ptr(vertex.position));
        }
    }
    glEnd();
}
