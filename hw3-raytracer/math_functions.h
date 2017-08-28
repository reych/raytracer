#ifndef MATH_FUNCTIONS_H
#define MATH_FUNCTIONS_H
#include "glm/glm.hpp"

typedef int sf_int;
const sf_int SPHERE_DOUBLE_ROOT = 2;
const sf_int SPHERE_UNEQUAL_ROOT = 1;
const sf_int SPHERE_NO_ROOT = 0;

float sphereIntersection(const glm::vec3 & ray0, const glm::vec3 & rayD, const glm::vec3 & spherePos, float radius, sf_int & sphereRoot);
float triangleIntersection(const glm::vec3 & ray0, const glm::vec3 & rayD, const glm::vec3 & posA,
    const glm::vec3 & posB, const glm::vec3 & posC, glm::vec3 & barycentricCoord);

glm::vec3 reflectionVector(const glm::vec3 & ray0, const glm::vec3 & rayD, const glm::vec3 & normal);

glm::vec3 calculate2DBarycentricCoord(const glm::vec3 & point, const glm::vec3 & posA, const glm::vec3 & posB, const glm::vec3 & posC);

float round2PlacesDecimal(float val);
float clamp(float value, float low, float high);

template <typename T>
glm::vec3 arrayToVec3(T arr[]) {
    return glm::vec3((float)arr[0], (float)arr[1], (float)arr[2]);
}

#endif
