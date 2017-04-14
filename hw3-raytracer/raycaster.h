#ifndef RAYCASTER_H
#define RAYCASTER_H
#include "glm/glm.hpp"

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

float calculateSphereIntersection(const glm::vec3 & ray0, const glm::vec3 & rayD, const glm::vec3 & spherePos, float radius);
float calculateTriangleIntersection(const glm::vec3 & ray0, const glm::vec3 & rayD, const glm::vec3 & posA, const glm::vec3 & posB, const glm::vec3 & posC);

glm::vec3 fireShadowRay(glm::vec3 pos, glm::vec3 viewer, glm::vec3 normal, float kd[], float ks[], float shininess,
    Light lights[], int numLights, Sphere spheres[], int numSpheres, Triangle triangles[], int numTriangles);

glm::vec3 calculate2DBarycentricCoord(const glm::vec3 & point, const glm::vec3 & posA, const glm::vec3 & posB, const glm::vec3 & posC);
float round2PlacesDecimal(float val);
float clamp(float value, float low, float high);

#endif
