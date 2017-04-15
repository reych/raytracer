#ifndef RAYTRACER_H
#define RAYTRACER_H
#include "glm/glm.hpp"

const int SPHERE = 0;
const int TRIANGLE = 1;
const int LIGHT = 2;
const int NOTHING = -1;

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

glm::vec3 localColor(const glm::vec3 & intersection, const glm::vec3 & viewer, int type, int index,
    Light lights[], int numLights, Sphere spheres[], int numSpheres, Triangle triangles[], int numTriangles, const glm::vec3 & barycentricCoord);
glm::vec3 fireShadowRay(glm::vec3 pos, glm::vec3 viewer, glm::vec3 normal, float kd[], float ks[], float shininess,
    Light lights[], int numLights, Sphere spheres[], int numSpheres, Triangle triangles[], int numTriangles);
glm::vec3 trace(const glm::vec3 & ray0, const glm::vec3 & pos, glm::vec3 viewer, int type, int index,
    Light lights[], int numLights, Sphere spheres[], int numSpheres, Triangle triangles[], int numTriangles, int recurseDepth);

glm::vec3 sphereNormal(const glm::vec3 & intersection, Sphere sphere);
glm::vec3 interpolateNormal(const glm::vec3 & barycentricCoord, Triangle triangle);

glm::vec3 closestIntersection(const glm::vec3 & ray0, const glm::vec3 & rayD,
    Sphere spheres[], int numSpheres, Triangle triangles[], int numTriangles, int & type, int & index, glm::vec3 & barycentricCoord);

#endif
