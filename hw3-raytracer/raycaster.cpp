#include "raycaster.h"
#include <stdlib.h>
#include <iostream>
using namespace std;

const int XY_PLANE = 0;
const int XZ_PLANE = 1;
const int YZ_PLANE = 2;

// Given ray origin 'ray0', normalized ray direction 'rayD', sphere center position 'spherePos', and sphere radius 'radius'.
// Return the minimum t value of the intersection.
float calculateSphereIntersection(const glm::vec3 & ray0, const glm::vec3 & rayD, const glm::vec3 & spherePos, float radius) {
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
        if(t0>0 && t1>0) {
            return fmin(t0, t1);
        }

        // Ignore negative t value (intersection is behind ray).
        if(t0<0.00001 && t1<0.00001) {
            return -1;
        }

        // One value is positive
        if(t0 < 0) {
            return t1;
        } else {
            return t0;
        }
    }
}

// Given ray origin 'ray0', ray direction 'rayD', triangle position 'posA, posB, posC', triangle normal 'normal'.
// Return the t parameter of the intersection point or -1 if it's not valid.
float calculateTriangleIntersection(const glm::vec3 & ray0, const glm::vec3 & rayD, const glm::vec3 & posA, const glm::vec3 & posB, const glm::vec3 & posC) {
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
    if(t <= 0.00015) {
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

    // Test alpha, beta, gamma values, return t if point within triangle.
    if((alpha >= 0 && alpha <= 1) && (beta >=0 && beta <=1) && (gamma >=0 && gamma <= 1)
        && alpha+beta+gamma > 0.999 && alpha+beta+gamma < 1.00001) {

        return t;
    }

    return -1;

}

// Fire shadow ray from point to light source.
// Return the color of the object RGB.
glm::vec3 fireShadowRay(glm::vec3 pos, glm::vec3 viewer, glm::vec3 normal, float kd[],float ks[], float shininess,
    Light lights[], int numLights, Sphere spheres[], int numSpheres, Triangle triangles[], int numTriangles) {
    glm::vec3 color = glm::vec3(0,0,0);
    // Calculate shadow ray for all lights.
    for(int i=0; i<numLights; i++) {
        bool blocked = false;
        // Create a ray.
        glm::vec3 ray0 = glm::vec3(pos.x, pos.y, pos.z);
        glm::vec3 l = glm::vec3(lights[i].position[0]-pos.x, lights[i].position[1]-pos.y, lights[i].position[2]-pos.z);
        l = glm::normalize(l); // Vector to light.

        // Check sphere intersections.
        for(int j=0; j<numSpheres; j++) {
            glm::vec3 spherePos = glm::vec3(spheres[j].position[0], spheres[j].position[1], spheres[j].position[2]);
            float radius = spheres[j].radius;
            float t = calculateSphereIntersection(ray0, l, spherePos, radius);
            if(t > 0.00015) {
                blocked = true;
            }
        }

        // Check triangle intersections.
        if(!blocked) {
            for(int j=0; j<numTriangles; j++) {
                glm::vec3 posA = glm::vec3(triangles[j].v[0].position[0], triangles[j].v[0].position[1], triangles[j].v[0].position[2]);
                glm::vec3 posB = glm::vec3(triangles[j].v[1].position[0], triangles[j].v[1].position[1], triangles[j].v[1].position[2]);
                glm::vec3 posC = glm::vec3(triangles[j].v[2].position[0], triangles[j].v[2].position[1], triangles[j].v[2].position[2]);
                float t = calculateTriangleIntersection(ray0, l, posA, posB, posC);
                if(t > 0.00015) {
                    blocked = true;
                }
            }
        }

        // If not blocked, calculate the Phong lighting.
        if(!blocked) {
            normal = glm::normalize(normal);
            viewer = glm::normalize(viewer); // Vector to viewer.
            float lDotN = glm::dot(l, normal);
            glm::vec3 r = glm::vec3(-l.x + 2*lDotN*normal.x, -l.y + 2*lDotN*normal.y, -l.z + 2*lDotN*normal.z); // Reflection vector.
            r = glm::normalize(r);

            float specularComp = clamp(glm::dot(r,viewer), 0, 1);
            float diffuseComp = clamp(glm::dot(l, normal), 0, 1);
            float red = lights[i].color[0] * ( kd[0]*(diffuseComp) + ks[0]*(pow(specularComp, shininess)) );
            float green = lights[i].color[1] * ( kd[1]*(diffuseComp) + ks[1]*(pow(specularComp, shininess)) );
            float blue = lights[i].color[2] * ( kd[2]*(diffuseComp) + ks[2]*(pow(specularComp, shininess)) );
            color = glm::vec3(color.x+red*255, color.y+green*255, color.z+blue*255);
        }

    }

    color = glm::vec3(clamp(color.x, 0, 255), clamp(color.y, 0, 255), clamp(color.z, 0, 255));
    return color;


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
