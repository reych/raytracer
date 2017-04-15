#include "math_functions.h"
#include <stdlib.h>
#include <iostream>
using namespace std;

// Given ray origin 'ray0', normalized ray direction 'rayD', sphere center position 'spherePos', and sphere radius 'radius'.
// Return the minimum t value of the intersection.
float sphereIntersection(const glm::vec3 & ray0, const glm::vec3 & rayD, const glm::vec3 & spherePos, float radius) {
    //rayD = glm::normalize(rayD);
    float b = 2*( rayD.x*(ray0.x-spherePos.x) + rayD.y*(ray0.y-spherePos.y) + rayD.z*(ray0.z-spherePos.z) );
    float c = (ray0.x-spherePos.x)*(ray0.x-spherePos.x) + (ray0.y-spherePos.y)*(ray0.y-spherePos.y) + (ray0.z-spherePos.z)*(ray0.z-spherePos.z) - radius*radius;

    float discriminant = b*b - 4*c;

    if(discriminant < 0.00001) {
        return -1;
    } else {
        float t0 = (-b + sqrt(discriminant))/2.0;
        float t1 = (-b - sqrt(discriminant))/2.0;

        // If neither are negative, return minimum of values.
        if(abs(t0)>1e-4 && abs(t1)>1e-4) {
            return fmin(t0, t1);
        }

        // Ignore negative t value (intersection is behind ray).
        if(t0<0.00001 && t1<0.00001) {
            return -1;
        }

        // One value is positive
        if(t0 < 0.00001) {
            return t1;
        } else {
            return t0;
        }
    }
}

// Given ray origin 'ray0', ray direction 'rayD', triangle position 'posA, posB, posC', triangle normal 'normal'.
// Return the t parameter of the intersection point or -1 if it's not valid. Output parameter barycentricCoord:
float triangleIntersection(const glm::vec3 & ray0, const glm::vec3 & rayD, const glm::vec3 & posA, const glm::vec3 & posB, const glm::vec3 & posC, glm::vec3 & barycentricCoord) {
    // Calculate normal.
    glm::vec3 BA = glm::vec3(posB.x-posA.x, posB.y-posA.y, posB.z-posA.z);
    glm::vec3 CA = glm::vec3(posC.x-posA.x, posC.y-posA.y, posC.z-posA.z);
    glm::vec3 normal = glm::cross(BA,CA);

    float normalDotRayD = glm::dot(normal, rayD);

    // Ray parallel to plane.
    if(abs(normalDotRayD) < 1e-6) {
        return -1;
    }

    // Otherwise, ray intersects plane.
    float d = -(glm::dot(normal, posA));
    float t = -(glm::dot(normal, ray0) + d) / normalDotRayD;

    // Intersection behind ray origin, exit.
    if(t <= 0.0001) {
        return -1;
    }

    // Otherwise, test if point is inside triangle.
    // --- [ Project triangle onto x=0, y=0, or z=0 plane. ] ---
    glm::vec3 xZero = glm::vec3(0, 1, 1); //yz plane
    glm::vec3 yZero = glm::vec3(1, 0, 1); //xz plane
    glm::vec3 zZero = glm::vec3(1, 1, 0); //xy plane

    // Take the angle between normal and plane.
    float xZeroAngle = abs(glm::dot(normal, xZero)); // yz plane.
    float yZeroAngle = abs(glm::dot(normal, yZero)); // xz plane.
    float zZeroAngle = abs(glm::dot(normal, zZero)); // xy plane.

    float minAngle = min(xZeroAngle, min(yZeroAngle, zZeroAngle)); //cout<<"min angle: "<<minAngle<<endl;
    glm::vec2 projectedA;
    glm::vec2 projectedB;
    glm::vec2 projectedC;
    glm::vec2 p;
    int projectionPlane;

    if(minAngle == xZeroAngle) {
        projectedA = glm::vec2(posA.y, posA.z);
        projectedB = glm::vec2(posB.y, posB.z);
        projectedC = glm::vec2(posC.y, posC.z);
        p = glm::vec2(ray0.y + rayD.y*t, ray0.z + rayD.z*t);
    } else if(minAngle == yZeroAngle) {
        projectedA = glm::vec2(posA.x, posA.z);
        projectedB = glm::vec2(posB.x, posB.z);
        projectedC = glm::vec2(posC.x, posC.z);
        p = glm::vec2(ray0.x + rayD.x*t, ray0.z + rayD.z*t);
    } else {
        projectedA = glm::vec2(posA.x, posA.y);
        projectedB = glm::vec2(posB.x, posB.y);
        projectedC = glm::vec2(posC.x, posC.y);
        p = glm::vec2(ray0.x + rayD.x*t, ray0.y + rayD.y*t);
    }

    // --- [ 2D test by Barycentric coordinates. ] ---
    // Get areas of inner triangles using projected points.
    float areaABC = .5 * ((projectedB.x - projectedA.x)*(projectedC.y - projectedA.y)
                            - (projectedC.x - projectedA.x)*(projectedB.y - projectedA.y));
    float areaApC = .5 * ((p.x - projectedA.x)*(projectedC.y - projectedA.y)
                            - (projectedC.x - projectedA.x)*(p.y - projectedA.y));
    float areaABp = .5 * ((projectedB.x - projectedA.x)*(p.y - projectedA.y)
                            - (p.x - projectedA.x)*(projectedB.y - projectedA.y));
    float areapBC = .5 * ((projectedB.x - p.x)*(projectedC.y - p.y)
                            - (projectedC.x - p.x)*(projectedB.y - p.y));

    // alpha, gamma, beta.
    float beta = areaApC/areaABC;
    float gamma = areaABp/areaABC;
    float alpha = areapBC/areaABC;

    barycentricCoord = glm::vec3(alpha, beta, gamma);

    // Test alpha, beta, gamma values, return t if point within triangle.
    if((alpha >= 0 && alpha <= 1) && (beta >=0 && beta <=1) && (gamma >=0 && gamma <= 1)
        && alpha+beta+gamma > 0.999 && alpha+beta+gamma < 1.00001) {

        return t;
    }

    return -1;

}

glm::vec3 reflectionVector(const glm::vec3 & ray0, const glm::vec3 & rayD, const glm::vec3 & normal) {
    glm::vec3 l = glm::vec3(rayD.x-ray0.x, rayD.y-ray0.y, rayD.z-ray0.z); // Incoming ray to intersection.
    float lDotN = glm::dot(l, normal);
    glm::vec3 r = glm::vec3(-l.x + 2*lDotN*normal.x, -l.y + 2*lDotN*normal.y, -l.z + 2*lDotN*normal.z); // Reflection vector.
    r = glm::normalize(r);
    return r;
}

// Given point and triangle vertices posA, posB, posC.
// Return Barycentric Coordinates alpha, beta, gamma.
glm::vec3 calculate2DBarycentricCoord(const glm::vec3 & point, const glm::vec3 & posA, const glm::vec3 & posB, const glm::vec3 & posC) {
    // Calculate normal.
    glm::vec3 BA = glm::vec3(posB.x-posA.x, posB.y-posA.y, posB.z-posA.z);
    glm::vec3 CA = glm::vec3(posC.x-posA.x, posC.y-posA.y, posC.z-posA.z);
    glm::vec3 normal = glm::cross(BA,CA);
    normal = glm::normalize(normal);

    // --- [ Project triangle onto x=0, y=0, or z=0 plane. ] ---
    glm::vec3 xZero = glm::vec3(0, 1, 1); //yz plane
    glm::vec3 yZero = glm::vec3(1, 0, 1); //xz plane
    glm::vec3 zZero = glm::vec3(1, 1, 0); //xy plane

    // Take the angle between normal and plane.
    float xZeroAngle = glm::dot(normal, xZero); // yz plane.
    float yZeroAngle = glm::dot(normal, yZero); // xz plane.
    float zZeroAngle = glm::dot(normal, zZero); // xy plane.

    float minAngle = min(xZeroAngle, min(yZeroAngle, zZeroAngle));
    glm::vec2 projectedA;
    glm::vec2 projectedB;
    glm::vec2 projectedC;
    glm::vec2 p;

    if(minAngle == xZeroAngle) {
        projectedA = glm::vec2(posA.y, posA.z);
        projectedB = glm::vec2(posB.y, posB.z);
        projectedC = glm::vec2(posC.y, posC.z);
        p = glm::vec2(point.y, point.z);
    } else if(minAngle == yZeroAngle) {
        projectedA = glm::vec2(posA.x, posA.z);
        projectedB = glm::vec2(posB.x, posB.z);
        projectedC = glm::vec2(posC.x, posC.z);
        p = glm::vec2(point.x, point.z);
    } else {
        projectedA = glm::vec2(posA.x, posA.y);
        projectedB = glm::vec2(posB.x, posB.y);
        projectedC = glm::vec2(posC.x, posC.y);
        p = glm::vec2(point.x, point.y);
    }

    // --- [ 2D test by Barycentric coordinates. ] ---
    // Get areas of inner triangles using projected points.
    float areaABC = .5 * ((projectedB.x - projectedA.x)*(projectedC.y - projectedA.y)
                            - (projectedC.x - projectedA.x)*(projectedB.y - projectedA.y));
    float areaApC = .5 * ((p.x - projectedA.x)*(projectedC.y - projectedA.y)
                            - (projectedC.x - projectedA.x)*(p.y - projectedA.y));
    float areaABp = .5 * ((projectedB.x - projectedA.x)*(p.y - projectedA.y)
                            - (p.x - projectedA.x)*(projectedB.y - projectedA.y));
    float areapBC = .5 * ((projectedB.x - p.x)*(projectedC.y - p.y)
                            - (projectedC.x - p.x)*(projectedB.y - p.y));

    // alpha, gamma, beta.
    float beta = areaApC/areaABC;
    float gamma = areaABp/areaABC;
    float alpha = areapBC/areaABC;

    glm::vec3 bCoord = glm::vec3(alpha, beta, gamma);
    return bCoord;
}

float round2PlacesDecimal(float val) {
    float nearest = roundf(val * 100) / 100;
    return nearest;
}

float clamp(float value, float low, float high) {
    if(value < low) {
        return low;
    } else if(value > high) {
        return high;
    }
    return value;
}
