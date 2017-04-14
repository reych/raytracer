/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Rena Chen
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include "glm/glm.hpp"
#include <iostream>
#include "raycaster.h"
using namespace std;

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

//MODIFY THIS FUNCTION
void draw_scene()
{
    // For every pixel, fire a ray to find first intersection (These are the visible vertices).
    // -->fire shadow ray to determine color.
    // Then put the color in the buffer.
    float pointX = -1.0*(1.0*WIDTH/HEIGHT*tan((fov*3.14159/180.0)/2.0));
    for(unsigned int x=0; x<WIDTH; x++) {
        glPointSize(2.0);
        glBegin(GL_POINTS);

        float pointY = -1.0* (tan((fov*3.14159/180.0)/2.0) );

        for(unsigned int y=0; y<HEIGHT; y++) {
            glm::vec3 ray0 = glm::vec3(0,0,0);
            glm::vec3 rayD = glm::normalize(glm::vec3(pointX, pointY, -1));
            int s_i = -1; // Sphere index.
            int t_i = -1; // Triangle index.
            float t_s = -1; // Sphere t parameter.
            float t_tr = -1; // Triangle t parameter.

            //Find intersections with spheres.
            for(int s=0; s<num_spheres; s++) {
                glm::vec3 pos = glm::vec3(spheres[s].position[0], spheres[s].position[1], spheres[s].position[2]);
                float radius = spheres[s].radius;
                float t_temp = calculateSphereIntersection(ray0, rayD, pos, radius);

                // Found intersection.
                if(t_temp > -1) {
                    // Get the closest point.
                    if(t_s > t_temp || t_s == -1) {
                        t_s = t_temp;
                        s_i = s;
                    }
                }
            }

            //Find intersections with triangles.
            for(int t=0; t<num_triangles; t++) {
                glm::vec3 posA = glm::vec3((float)triangles[t].v[0].position[0], (float)triangles[t].v[0].position[1], (float)triangles[t].v[0].position[2]);
                glm::vec3 posB = glm::vec3((float)triangles[t].v[1].position[0], (float)triangles[t].v[1].position[1], (float)triangles[t].v[1].position[2]);
                glm::vec3 posC = glm::vec3((float)triangles[t].v[2].position[0], (float)triangles[t].v[2].position[1], (float)triangles[t].v[2].position[2]);
                float t_temp = calculateTriangleIntersection(ray0, rayD, posA, posB, posC);

                // Found intersection.
                if(t_temp > -1) {
                    // Get the closest point.
                    if(t_tr > t_temp || t_tr == -1) {
                        t_tr = t_temp;
                        t_i = t;
                    }

                }
            }

            // Determine the closest intersection.
            bool doSphere = ((t_s != -1 && t_tr != -1) && (t_s < t_tr)) || (t_tr == -1 && t_s != -1);
            bool doTriangle = ((t_s != -1 && t_tr != -1) && (t_tr < t_s)) || (t_s == -1 && t_tr != -1);
            glm::vec3 color = glm::vec3(0, 0, 0); // Adjust ambient color here.
            if(doSphere) {
                // Get sphere color.
                glm::vec3 intersection = glm::vec3(ray0.x + rayD.x*t_s, ray0.y + rayD.y*t_s, ray0.z + rayD.z*t_s);
                glm::vec3 normal = (float)(1/spheres[s_i].radius) * glm::vec3(intersection.x-(float)spheres[s_i].position[0], intersection.y-(float)spheres[s_i].position[1], intersection.z-(float)spheres[s_i].position[2]);
                float ks[3] = {spheres[s_i].color_specular[0], spheres[s_i].color_specular[1], spheres[s_i].color_specular[2]};
                float kd[3] = {spheres[s_i].color_diffuse[0], spheres[s_i].color_diffuse[1], spheres[s_i].color_diffuse[2]};
                // TODO: negate if inside of sphere.

                glm::vec3 color_s = fireShadowRay(intersection, -rayD, normal, kd, ks,
                                                (float)spheres[s_i].shininess, lights, num_lights, spheres, num_spheres, triangles, num_triangles);
                color = glm::vec3(color.x+color_s.x, color.y+color_s.y, color.z+color_s.z);

            } else if(doTriangle){
                // Get triangle color.
                glm::vec3 intersection = glm::vec3(ray0.x + rayD.x*t_tr, ray0.y + rayD.y*t_tr, ray0.z + rayD.z*t_tr);
                // Interpolate normals, diffuse, normal, shininess.
                glm::vec3 projectedA = glm::vec3(triangles[t_i].v[0].position[0], triangles[t_i].v[0].position[1], triangles[t_i].v[0].position[2]);
                glm::vec3 projectedB = glm::vec3(triangles[t_i].v[1].position[0], triangles[t_i].v[1].position[1], triangles[t_i].v[1].position[2]);
                glm::vec3 projectedC = glm::vec3(triangles[t_i].v[2].position[0], triangles[t_i].v[2].position[1], triangles[t_i].v[2].position[2]);
                glm::vec3 bCoord = calculate2DBarycentricCoord(intersection, projectedA, projectedB, projectedC);
                float alpha = bCoord.x;
                float beta = bCoord.y;
                float gamma = bCoord.z;

                // triangle ABC,
                // normals nA, nB, nC,
                // alpha*nA, beta*nB, gamma*nC,
                // = alpha(nA.x, nA.y, nA.z), beta(nB.x, nB.y, nB.z), gamma*(nC.x, nC.y, nC.z);
                // normal = (alpha*nA.x + beta*nB.x + gamma*nC.x, alpha*nA.y + beta*nB.y + gamma*nC.y, alpha*nA.z + beta*nB.z + gamma*nC.z);
                double* nA = triangles[t_i].v[0].normal;
                double* nB = triangles[t_i].v[1].normal;
                double* nC = triangles[t_i].v[2].normal;
                glm::vec3 normal = glm::vec3((float)(alpha*nA[0] + beta*nB[0] + gamma*nC[0]),
                                                (float)(alpha*nA[1] + beta*nB[1] + gamma*nC[1]),
                                                (float)(alpha*nA[2] + beta*nB[2] + gamma*nC[2]));
                // diffuse = (alpha*dA.R + beta*dB.R + gamma*dC.R, alpha*dA.G + beta*dB.G + gamma*dC.G, alpha*dA.B + beta*dB.B + gamma*dC.B)
                double* dA = triangles[t_i].v[0].color_diffuse;
                double* dB = triangles[t_i].v[1].color_diffuse;
                double* dC = triangles[t_i].v[2].color_diffuse;
                float diffuse[3] = {(float)(alpha*dA[0] + beta*dB[0] + gamma*dC[0]),
                                        (float)(alpha*dA[1] + beta*dB[1] + gamma*dC[1]),
                                        (float)(alpha*dA[2] + beta*dB[2] + gamma*dC[2])};
                double* spA = triangles[t_i].v[0].color_specular;
                double* spB = triangles[t_i].v[1].color_specular;
                double* spC = triangles[t_i].v[2].color_specular;
                float specular[3] = {(float)(alpha*spA[0] + beta*spB[0] + gamma*spC[0]),
                                        (float)(alpha*spA[1] + beta*spB[1] + gamma*spC[1]),
                                        (float)(alpha*spA[2] + beta*spB[2] + gamma*spC[2])};
                double shA = triangles[t_i].v[0].shininess;
                double shB = triangles[t_i].v[1].shininess;
                double shC = triangles[t_i].v[2].shininess;
                float shininess = (float)(alpha*shA + beta*shB + gamma*shC);

                glm::vec3 color_tr = fireShadowRay(intersection, -rayD, normal, diffuse, specular, shininess,
                                                lights, num_lights, spheres, num_spheres, triangles, num_triangles);
                color = glm::vec3(color.x+color_tr.x, color.y+color_tr.y, color.z+color_tr.z);

            } else {
                color = glm::vec3(255, 255, 255);
            }

            plot_pixel(x, y, clamp(color.x,0,255), clamp(color.y,0,255), clamp(color.z,0,255));
            pointY += (2.0*tan((fov*3.14159/180.0)/2.0)) /  HEIGHT;
        }
        glEnd();
        glFlush();
        pointX += (2.0*WIDTH/HEIGHT*tan((fov*3.14159/180.0)/2.0)) / WIDTH;
    }

  printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
