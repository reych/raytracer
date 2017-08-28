#include "raytracer.h"
#include "math_functions.h"
#include <stdlib.h>
#include <iostream>
using namespace std;

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
            sf_int sphereRootFlag;
            float t = sphereIntersection(ray0, l, spherePos, radius, sphereRootFlag);
            if(t > 0.00015 && sphereRootFlag != SPHERE_DOUBLE_ROOT) {
                blocked = true;
                break;
            }
        }
        // Check triangle intersections.
        if(!blocked) {
            for(int j=0; j<numTriangles; j++) {
                glm::vec3 posA = glm::vec3(triangles[j].v[0].position[0], triangles[j].v[0].position[1], triangles[j].v[0].position[2]);
                glm::vec3 posB = glm::vec3(triangles[j].v[1].position[0], triangles[j].v[1].position[1], triangles[j].v[1].position[2]);
                glm::vec3 posC = glm::vec3(triangles[j].v[2].position[0], triangles[j].v[2].position[1], triangles[j].v[2].position[2]);
                glm::vec3 barycentricCoord;
                float t = triangleIntersection(ray0, l, posA, posB, posC, barycentricCoord);
                if(t > 0.00015) {
                    blocked = true;
                    break;
                }
            }
        }

        // If not blocked, calculate the Phong lighting.
        if(!blocked) {
            normal = glm::normalize(normal);
            viewer = -pos;
            viewer = glm::normalize(viewer); // Vector to viewer.
            float lDotN = glm::dot(l, normal);
            glm::vec3 r = glm::vec3(-l.x + 2*lDotN*normal.x, -l.y + 2*lDotN*normal.y, -l.z + 2*lDotN*normal.z); // Reflection vector.
            //glm::vec3 r = -glm::reflect(l, normal);
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

// Return color from reflection and transmission rays.
glm::vec3 trace(const glm::vec3 & ray0, const glm::vec3 & rayD, glm::vec3 camera, int type, int index,
    Light lights[], int numLights, Sphere spheres[], int numSpheres, Triangle triangles[], int numTriangles, int recurseDepth) {

    glm::vec3 local, reflected, transmitted;

    // If reached end of recursion.
    if (recurseDepth == 0)
        return glm::vec3(0, 0, 0);

    // Find intersection.
    int type_intersect = -1;
    int index_intersect = -1;
    glm::vec3 barycentricCoord; // Output parameter for triangle interpolation functions.
    glm::vec3 intersection = closestIntersection(ray0, rayD, spheres, numSpheres, triangles, numTriangles, type_intersect, index_intersect, barycentricCoord);

    // If hit light, return light color.
    if(type_intersect == LIGHT)
        return (arrayToVec3(lights[index_intersect].color));
    // If no intersection, return black.
    if(type_intersect == NOTHING)
        return glm::vec3(0,0,0);

    glm::vec3 normal;
    if(type_intersect == SPHERE)
        normal = sphereNormal(rayD, spheres[index_intersect]);
    else if(type_intersect ==TRIANGLE)
        normal = interpolateNormal(barycentricCoord, triangles[index_intersect]);

    glm::vec3 l = glm::normalize(glm::vec3(ray0.x-rayD.x, ray0.y-rayD.y, ray0.z-rayD.z));

    normal = glm::normalize(normal); //Should i do this?
    float lDotN = glm::dot(l, normal);

    glm::vec3 r = glm::vec3(-l.x + 2*lDotN*normal.x, -l.y + 2*lDotN*normal.y, -l.z + 2*lDotN*normal.z); // Reflection vector.
    r = glm::normalize(r);
    glm::vec3 viewer = glm::vec3(camera.x-intersection.x, camera.y-intersection.y, camera.z-intersection.z);

    int upperType = NOTHING;
    int lowerType = NOTHING;
    // If same type, and same index, then same object, no refraction because going through same medium
    if(type == type_intersect) {
        if(index == index_intersect) {
            upperType = type;
            lowerType = type;
        } else {
            upperType = NOTHING;
            lowerType = type_intersect;
        }

    }
    // if not same object, then assume there's air in between, so first material is air.
    else {
        upperType = NOTHING;
        lowerType = type_intersect;
    }
    glm::vec3 t = glm::normalize(transmit(normal, l, upperType, lowerType)); // Get transmitted angle.

    local = localColor(intersection, viewer, type_intersect, index_intersect, lights, numLights, spheres, numSpheres, triangles, numTriangles, barycentricCoord);
    reflected = trace(intersection, r, camera, type_intersect, index_intersect, lights, numLights, spheres, numSpheres, triangles, numTriangles, recurseDepth - 1);
    //transmitted = trace(intersection, t, camera, type_intersect, index_intersect, lights, numLights, spheres, numSpheres, triangles, numTriangles, recurseDepth - 1);

    float discount = .3;
    float t_discount = .3;
    float alpha = 0.5;
    return glm::vec3(local.x + discount*reflected.x, local.y + discount*reflected.y, local.z + discount*reflected.z);
    //return glm::vec3(local.x + discount*reflected.x + t_discount*transmitted.x, local.y + discount*reflected.y + t_discount*transmitted.y, local.z + discount*reflected.z + t_discount*transmitted.z);
    //return glm::vec3(local.x + t_discount*transmitted.x, local.y + t_discount*transmitted.y, local.z + t_discount*transmitted.z);


}

// Given normal, incoming angle l, upper material 'upperType', and lower material 'lowerMaterial'.
// Return vector direction of transmission.
glm::vec3 transmit(const glm::vec3 & normal, const glm::vec3 & l, int upperType, int lowerType) {
    // Set some refractive index values.
    float nl = 1.0;// index of upper material.
    float nt = 1.0;// index of lower material.

    if (upperType == SPHERE) {
        nl = 1.33; // Water.
    } else if (upperType == TRIANGLE) {
        nl = 1.49;
    }
    if (lowerType == SPHERE) {
        nt = 1.33; // Water.
    } else if (upperType == TRIANGLE) {
        nt = 1.49;
    }
    float n = nt/nl;

    return refract(normal, l, n);

}

// Given 'intersection' point, 'viewer' vector, 'type' of object, 'index' of object, and arrays of lights and objects.
// Return local color. (Wrapper around fireShadowRay.)
glm::vec3 localColor(const glm::vec3 & intersection, const glm::vec3 & viewer, int type, int index,
    Light lights[], int numLights, Sphere spheres[], int numSpheres, Triangle triangles[], int numTriangles, const glm::vec3 & barycentricCoord) {

    if(type == SPHERE) {
        // Get sphere color.
        glm::vec3 normal = (float)(1/spheres[index].radius) * glm::vec3(intersection.x-(float)spheres[index].position[0], intersection.y-(float)spheres[index].position[1], intersection.z-(float)spheres[index].position[2]);
        float ks[3] = {spheres[index].color_specular[0], spheres[index].color_specular[1], spheres[index].color_specular[2]};
        float kd[3] = {spheres[index].color_diffuse[0], spheres[index].color_diffuse[1], spheres[index].color_diffuse[2]};
        // TODO: negate if inside of sphere.

        glm::vec3 color_s = fireShadowRay(intersection, viewer, normal, kd, ks,
            (float)spheres[index].shininess, lights, numLights, spheres, numSpheres, triangles, numTriangles);
        return color_s;
    } else if(type == TRIANGLE) {
        // Get triangle color.
        // Interpolate normals, diffuse, normal, shininess.
        glm::vec3 projectedA = glm::vec3(triangles[index].v[0].position[0], triangles[index].v[0].position[1], triangles[index].v[0].position[2]);
        glm::vec3 projectedB = glm::vec3(triangles[index].v[1].position[0], triangles[index].v[1].position[1], triangles[index].v[1].position[2]);
        glm::vec3 projectedC = glm::vec3(triangles[index].v[2].position[0], triangles[index].v[2].position[1], triangles[index].v[2].position[2]);
        //glm::vec3 bCoord = calculate2DBarycentricCoord(intersection, projectedA, projectedB, projectedC);
        float alpha = barycentricCoord.x;
        float beta = barycentricCoord.y;
        float gamma = barycentricCoord.z;

        // triangle ABC,
        // normals nA, nB, nC,
        // alpha*nA, beta*nB, gamma*nC,
        // = alpha(nA.x, nA.y, nA.z), beta(nB.x, nB.y, nB.z), gamma*(nC.x, nC.y, nC.z);
        // normal = (alpha*nA.x + beta*nB.x + gamma*nC.x, alpha*nA.y + beta*nB.y + gamma*nC.y, alpha*nA.z + beta*nB.z + gamma*nC.z);
        double* nA = triangles[index].v[0].normal;
        double* nB = triangles[index].v[1].normal;
        double* nC = triangles[index].v[2].normal;
        glm::vec3 normal = glm::vec3((float)(alpha*nA[0] + beta*nB[0] + gamma*nC[0]),
                                        (float)(alpha*nA[1] + beta*nB[1] + gamma*nC[1]),
                                        (float)(alpha*nA[2] + beta*nB[2] + gamma*nC[2]));
        // diffuse = (alpha*dA.R + beta*dB.R + gamma*dC.R, alpha*dA.G + beta*dB.G + gamma*dC.G, alpha*dA.B + beta*dB.B + gamma*dC.B)
        double* dA = triangles[index].v[0].color_diffuse;
        double* dB = triangles[index].v[1].color_diffuse;
        double* dC = triangles[index].v[2].color_diffuse;
        float diffuse[3] = {(float)(alpha*dA[0] + beta*dB[0] + gamma*dC[0]),
                                (float)(alpha*dA[1] + beta*dB[1] + gamma*dC[1]),
                                (float)(alpha*dA[2] + beta*dB[2] + gamma*dC[2])};
        double* spA = triangles[index].v[0].color_specular;
        double* spB = triangles[index].v[1].color_specular;
        double* spC = triangles[index].v[2].color_specular;
        float specular[3] = {(float)(alpha*spA[0] + beta*spB[0] + gamma*spC[0]),
                                (float)(alpha*spA[1] + beta*spB[1] + gamma*spC[1]),
                                (float)(alpha*spA[2] + beta*spB[2] + gamma*spC[2])};
        double shA = triangles[index].v[0].shininess;
        double shB = triangles[index].v[1].shininess;
        double shC = triangles[index].v[2].shininess;
        float shininess = (float)(alpha*shA + beta*shB + gamma*shC);

        glm::vec3 color_tr = fireShadowRay(intersection, viewer, normal, diffuse, specular, shininess,
                                        lights, numLights, spheres, numSpheres, triangles, numTriangles);
        return color_tr;
    }
    return glm::vec3(0,0,0);
}

glm::vec3 sphereNormal(const glm::vec3 & intersection, Sphere sphere) {
    glm::vec3 normal = (float)(1.0/sphere.radius) * glm::vec3(intersection.x-(float)sphere.position[0], intersection.y-(float)sphere.position[1], intersection.z-(float)sphere.position[2]);
    return normal;
}
glm::vec3 interpolateNormal(const glm::vec3 & barycentricCoord, Triangle triangle) {
    glm::vec3 projectedA = glm::vec3(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
    glm::vec3 projectedB = glm::vec3(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
    glm::vec3 projectedC = glm::vec3(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);
    float alpha = barycentricCoord.x;
    float beta = barycentricCoord.y;
    float gamma = barycentricCoord.z;

    double* nA = triangle.v[0].normal;
    double* nB = triangle.v[1].normal;
    double* nC = triangle.v[2].normal;
    glm::vec3 normal = glm::vec3((float)(alpha*nA[0] + beta*nB[0] + gamma*nC[0]),
                                    (float)(alpha*nA[1] + beta*nB[1] + gamma*nC[1]),
                                    (float)(alpha*nA[2] + beta*nB[2] + gamma*nC[2]));
    return normal;
}

// Calculates the closest intersection relative to vector ray0 -> rayD.
// Output parameters: intersection parameter 't', type of object 'type', index of object 'index'.
// Return intersection point.
glm::vec3 closestIntersection(const glm::vec3 & ray0, const glm::vec3 & rayD,
    Sphere spheres[], int numSpheres, Triangle triangles[], int numTriangles, int & type, int & index, glm::vec3 & barycentricCoord) {
        int s_i = -1; // Sphere index.
        int t_i = -1; // Triangle index.
        float t_s = -1; // Sphere t parameter.
        float t_tr = -1; // Triangle t parameter.
        float t;

        //Find intersections with spheres.
        for(int s=0; s<numSpheres; s++) {
            glm::vec3 pos = glm::vec3(spheres[s].position[0], spheres[s].position[1], spheres[s].position[2]);
            float radius = spheres[s].radius;
            sf_int sphereRootFlag;
            float t_temp = sphereIntersection(ray0, rayD, pos, radius, sphereRootFlag); // TODO: May be intersecting with self...need to get closest intersection that is not a doubleroot.

            // Found intersection.
            if(t_temp > -1 && sphereRootFlag != SPHERE_DOUBLE_ROOT) {
                // Get the closest point.
                if(t_s > t_temp || t_s == -1) {
                    t_s = t_temp;
                    s_i = s;
                }
            }
        }
        //Find intersections with triangles.
        for(int t=0; t<numTriangles; t++) {
            glm::vec3 posA = glm::vec3((float)triangles[t].v[0].position[0], (float)triangles[t].v[0].position[1], (float)triangles[t].v[0].position[2]);
            glm::vec3 posB = glm::vec3((float)triangles[t].v[1].position[0], (float)triangles[t].v[1].position[1], (float)triangles[t].v[1].position[2]);
            glm::vec3 posC = glm::vec3((float)triangles[t].v[2].position[0], (float)triangles[t].v[2].position[1], (float)triangles[t].v[2].position[2]);
            glm::vec3 bCoord;
            float t_temp = triangleIntersection(ray0, rayD, posA, posB, posC, bCoord);

            // Found intersection.
            if(t_temp > -1) {
                // Get the closest point.
                if(t_tr > t_temp || t_tr == -1) {
                    t_tr = t_temp;
                    t_i = t;
                    barycentricCoord = bCoord; // Only update barycenric coordinates for relevant triangle.
                }

            }
        }

        // Determine the closest intersection.
        bool doSphere = ((t_s != -1 && t_tr != -1) && (t_s < t_tr)) || (t_tr == -1 && t_s != -1);
        bool doTriangle = ((t_s != -1 && t_tr != -1) && (t_tr < t_s)) || (t_s == -1 && t_tr != -1);
        if(doSphere) {
            type = SPHERE;
            t = t_s;
            index = s_i;
            barycentricCoord = glm::vec3(1,1,1);
        } else if(doTriangle) {
            type = TRIANGLE;
            t = t_tr;
            index = t_i;
        } else {
            type = NOTHING;
            t = -1;
            index = -1;
            barycentricCoord = glm::vec3(0,0,0);
        }
        glm::vec3 intersection = glm::vec3(ray0.x + rayD.x*t, ray0.y + rayD.y*t, ray0.z + rayD.z*t);
        return intersection;
}
