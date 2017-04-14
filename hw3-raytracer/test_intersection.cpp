#include "raycaster.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
using namespace std;

int main() {
    // --------- [ Test circle intersection ] ---------
    // Circle intersects:
    glm::vec3 ray0 = glm::vec3(0,0,0);
    glm::vec3 rayD = glm::vec3(0, 4, -7);
    rayD = glm::normalize(rayD);
    glm::vec3 spherePos = glm::vec3(0, 4, -7);
    float radius = 2;
    float answerX = 0; float answerY = 3.01; float answerZ = -5.26;
    float t = calculateSphereIntersection(ray0, rayD, spherePos, radius);

    float actualX = round2PlacesDecimal(ray0.x + rayD.x*t);
    float actualY = round2PlacesDecimal(ray0.y + rayD.y*t);
    float actualZ = round2PlacesDecimal(ray0.z + rayD.z*t);

    if(actualX == answerX && actualY == answerY && actualZ == answerZ) {
        cout<<"PASS sphere intersects "<<endl;
    } else {
        cout<<"FAIL sphere intersects "<<endl;
        cout<<"    expected: "<<answerX<<", "<<answerY<<", "<<answerZ<<endl;
    }
    cout<<"    Actual answer: "<<actualX<<", "<<actualY<<", "<<actualZ<<endl;
    cout<<"    t: "<<t<<endl;

    // Ray misses sphere:
    rayD = glm::vec3(0, 10, -7);
    rayD = glm::normalize(rayD);
    t = calculateSphereIntersection(ray0, rayD, spherePos, radius);

    if(t == -1) {
        cout<<"PASS sphere miss "<<endl;
    } else {
        cout<<"FAIL sphere miss "<<endl;
        cout<<"    expected: "<<-1<<endl;
    }
    cout<<"    Actual answer: "<<t<<endl;

    // Sphere behind ray origin:
    ray0 = glm::vec3(0,0,0);
    rayD = glm::vec3(0, 4, -7);
    spherePos = glm::vec3(0, 4, 7);

    t = calculateSphereIntersection(ray0, rayD, spherePos, radius);

    if(t == -1) {
        cout<<"PASS sphere behind ray "<<endl;
    } else {
        cout<<"FAIL sphere behind ray "<<endl;
        cout<<"    expected: "<<-1<<endl;
    }
    cout<<"    Actual answer: "<<t<<endl;

    // Ray origin in sphere:
    ray0 = glm::vec3(0,0,0);
    rayD = glm::vec3(0, 4, -7);
    rayD = glm::normalize(rayD);
    spherePos = glm::vec3(0, 0, 0);
    answerX = 0; answerY = 0.99; answerZ = -1.74;

    t = calculateSphereIntersection(ray0, rayD, spherePos, radius);

    actualX = round2PlacesDecimal(ray0.x + rayD.x*t);
    actualY = round2PlacesDecimal(ray0.y + rayD.y*t);
    actualZ = round2PlacesDecimal(ray0.z + rayD.z*t);

    if(actualX == answerX && actualY == answerY && actualZ == answerZ) {
        cout<<"PASS ray origin in sphere "<<endl;
    } else {
        cout<<"FAIL ray origin in sphere "<<endl;
        cout<<"    expected: "<<answerX<<", "<<answerY<<", "<<answerZ<<endl;
    }
    cout<<"    Actual answer: "<<actualX<<", "<<actualY<<", "<<actualZ<<endl;
    cout<<"    t: "<<t<<endl;

    // --------- [ Test triangle intersection ] ---------
    ray0 = glm::vec3(0,0,0);
    rayD = glm::vec3(0, 4, -7);
    rayD = glm::normalize(rayD);
    // Triangle intersects:
    glm::vec3 trianglePosA = glm::vec3(-4, 0, -7);
    glm::vec3 trianglePosB = glm::vec3(0, 7, -7);
    glm::vec3 trianglePosC = glm::vec3(4, 0, -7);
    answerX = 0; answerY = 4; answerZ = -7;

    t = calculateTriangleIntersection(ray0, rayD, trianglePosA, trianglePosB, trianglePosC);

    actualX = round2PlacesDecimal(ray0.x + rayD.x*t);
    actualY = round2PlacesDecimal(ray0.y + rayD.y*t);
    actualZ = round2PlacesDecimal(ray0.z + rayD.z*t);

    if(actualX == answerX && actualY == answerY && actualZ == answerZ) {
        cout<<"PASS triangle intersection 1 "<<endl;
    } else {
        cout<<"FAIL triangle intersection 1 "<<endl;
        cout<<"    expected: "<<answerX<<", "<<answerY<<", "<<answerZ<<endl;
    }
    cout<<"    Actual answer: "<<actualX<<", "<<actualY<<", "<<actualZ<<endl;
    cout<<"    t: "<<t<<endl;

    // Triangle intersects 2:
    trianglePosA = glm::vec3(-4, 0, -7);
    trianglePosB = glm::vec3(0, 7, -8);
    trianglePosC = glm::vec3(4, 0, -9);
    answerX = 0; answerY = 4.57; answerZ = -8;

    t = calculateTriangleIntersection(ray0, rayD, trianglePosA, trianglePosB, trianglePosC);

    actualX = round2PlacesDecimal(ray0.x + rayD.x*t);
    actualY = round2PlacesDecimal(ray0.y + rayD.y*t);
    actualZ = round2PlacesDecimal(ray0.z + rayD.z*t);

    if(actualX == answerX && actualY == answerY && actualZ == answerZ) {
        cout<<"PASS triangle intersection 2 "<<endl;
    } else {
        cout<<"FAIL triangle intersection 2 "<<endl;
        cout<<"    expected: "<<answerX<<", "<<answerY<<", "<<answerZ<<endl;
    }
    cout<<"    Actual answer: "<<actualX<<", "<<actualY<<", "<<actualZ<<endl;
    cout<<"    t: "<<t<<endl;

    // Triangle parallel:
    ray0 = glm::vec3(0,0,0);
    rayD = glm::vec3(0, 0, -7);
    rayD = glm::normalize(rayD);
    trianglePosA = glm::vec3(0, 4, -7);
    trianglePosB = glm::vec3(1, 4, -8);
    trianglePosC = glm::vec3(2, 4, -10);

    t = calculateTriangleIntersection(ray0, rayD, trianglePosA, trianglePosB, trianglePosC);

    if(t == -1) {
        cout<<"PASS triangle parallel "<<endl;
    } else {
        cout<<"FAIL triangle intersection 1 "<<endl;
        cout<<"    expected: "<<-1<<endl;
    }
    cout<<"    Actual answer: "<<t<<endl;

    // --------- [ Test Barycentric Coord ] ---------
    trianglePosA = glm::vec3(-4, 0, -7);
    trianglePosB = glm::vec3(0, 7, -8);
    trianglePosC = glm::vec3(4, 0, -9);
    glm::vec3 point = glm::vec3(0, 4.57, -8);
    glm::vec3 barycentricCoord = calculate2DBarycentricCoord(point, trianglePosA, trianglePosB, trianglePosC);
    cout<<"Barycentric coordinates "<<endl;
    cout<<"    "<<barycentricCoord.x<<", "<<barycentricCoord.y<<", "<<barycentricCoord.z<<endl;
    cout<<"    alpha+beta+gamma = "<<barycentricCoord.x + barycentricCoord.y + barycentricCoord.z<<endl;

    trianglePosA = glm::vec3(-4, 0, -7);
    trianglePosB = glm::vec3(0, 7, -7);
    trianglePosC = glm::vec3(4, 0, -7);
    point = glm::vec3(0, 4, -7);
    barycentricCoord = calculate2DBarycentricCoord(point, trianglePosA, trianglePosB, trianglePosC);
    cout<<"Barycentric coordinates "<<endl;
    cout<<"    "<<barycentricCoord.x<<", "<<barycentricCoord.y<<", "<<barycentricCoord.z<<endl;
    cout<<"    alpha+beta+gamma = "<<barycentricCoord.x + barycentricCoord.y + barycentricCoord.z<<endl;



}
