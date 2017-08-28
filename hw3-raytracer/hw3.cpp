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
#include "raytracer.h"
#include "math_functions.h"
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

#define ssFactor 1.0

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
    float pixelWidth = (2.0*WIDTH/HEIGHT*tan((fov*3.14159/180.0)/2.0)) / WIDTH;
    float pixelHeight = (2.0*tan((fov*3.14159/180.0)/2.0)) /  HEIGHT;

    glm::vec3 ray0 = glm::vec3(0,0,0);

    float pointX = -1.0*(1.0*WIDTH/HEIGHT*tan((fov*3.14159/180.0)/2.0));

    for(unsigned int x=0; x<WIDTH; x++) {
        glPointSize(2.0);
        glBegin(GL_POINTS);

        float pointY = -1.0* (tan((fov*3.14159/180.0)/2.0) );

        for(unsigned int y=0; y<HEIGHT; y++) {
            // Super sample
            glm::vec3 pixelColor = glm::vec3(0, 0, 0); // Adjust ambient color here.

            float subPointX = pointX + pixelWidth/(2 * ssFactor);
            float subPointY = pointY + pixelHeight/(2 * ssFactor);

            while(subPointX < pointX + pixelWidth) {
                while(subPointY < pointY + pixelHeight) {

                    glm::vec3 rayD = glm::normalize(glm::vec3(subPointX, subPointY, -1));
                    int type = NOTHING;
                    int index = -1;
                    glm::vec3 barycentricCoord;
                    glm::vec3 camera = glm::vec3(0,0,0); // camera position

                    glm::vec3 color = trace(ray0, rayD, camera, type, index,
                        lights, num_lights, spheres, num_spheres, triangles, num_triangles, 5);

                    pixelColor.x += color.x;
                    pixelColor.y += color.y;
                    pixelColor.z += color.z;

                    subPointY += pixelHeight/ssFactor;
                }
                subPointX += pixelWidth/ssFactor;
            }

            pixelColor.x /= ssFactor;
            pixelColor.y /= ssFactor;
            pixelColor.z /= ssFactor;

            plot_pixel(x, y, clamp(pixelColor.x,10,255), clamp(pixelColor.y,10,255), clamp(pixelColor.z,10,255));


            /* ----------- Original non multi sample ray tracing --------------- */
            // glm::vec3 ray0 = glm::vec3(0,0,0);
            // glm::vec3 rayD = glm::normalize(glm::vec3(pointX, pointY, -1));
            // int type = NOTHING;
            // int index = -1;
            // glm::vec3 barycentricCoord;
            // glm::vec3 camera = glm::vec3(0,0,0); // camera position
            // //glm::vec3 intersection = closestIntersection(ray0, rayD, spheres, num_spheres, triangles, num_triangles, type, index, barycentricCoord);
            // glm::vec3 color = glm::vec3(0, 0, 0); // Adjust ambient color here.
            // //glm::vec3 color_addition = glm::vec3(0, 0, 0);
            //
            // //if(type == SPHERE || type == TRIANGLE) {
            //     //color = localColor(intersection, -rayD, type, index, lights, num_lights, spheres, num_spheres, triangles, num_triangles, barycentricCoord);
            //     color = trace(ray0, rayD, camera, type, index,
            //         lights, num_lights, spheres, num_spheres, triangles, num_triangles, 20);
            //     //color = glm::vec3(0.5*color.x + 0.5*color_addition.x, 0.5*color.y+0.5*color_addition.y, 0.5*color.z + 0.5*color_addition.z);
            // //} else {
            //     //color = glm::vec3(255, 255, 255);
            // //}
            //
            // plot_pixel(x, y, clamp(color.x,10,255), clamp(color.y,10,255), clamp(color.z,10,255));
            // //plot_pixel(x, y, color_addition.x, color_addition.y, color_addition.z);
            // //plot_pixel(x, y, clamp(color_addition.x,0,255), clamp(color_addition.y,0,255), clamp(color_addition.z,0,255));
            // //pointY += (2.0*tan((fov*3.14159/180.0)/2.0)) /  HEIGHT;

            // Increment
            pointY += pixelHeight;
        }
        glEnd();
        glFlush();
        //pointX += (2.0*WIDTH/HEIGHT*tan((fov*3.14159/180.0)/2.0)) / WIDTH;
        pointX += pixelWidth;
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
