/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Nate Boxer
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
#include <cmath>
#include <cassert>
#include <limits>
#include <algorithm>

#include "point.h"
#include "utils.h"
#include "objects.h"

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

char *filename = NULL;

// The different display modes.
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

// While solving the homework, it is useful to make the below values smaller for debugging purposes.
// The still images that you need to submit with the homework should be at the below resolution (640x480).
// However, for your own purposes, after you have solved the homework, you can increase those values to obtain higher-resolution images.
#define WIDTH 640
#define HEIGHT 480

// The field of view of the camera, in degrees.
#define fov 60.0

// Buffer to store the image when saving it to a JPEG.
// FIX: needs to be height then width, not width then height
unsigned char buffer[HEIGHT][WIDTH][3];

// buffer to hold image as we do ray tracing
float image[WIDTH][HEIGHT][3];

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

Point camera(0, 0, 0);
Point background = Point(255, 255, 255).convert8BitToFloat();

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

unsigned long counter = 0;

typedef enum ObjectType
{
  SPHERE,
  TRIANGLE,
  LIGHT,
  NONE
} ObjectType;

struct Intersection
{
  double t;
  int i;
  Point b_coords;
  ObjectType type;

  Intersection(double t, int i, Point b, ObjectType type) : t(t), i(i), b_coords(b), type(type) {}

  bool operator==(const Intersection &rhs)
  {
    return t == rhs.t && i == rhs.i && b_coords == rhs.b_coords && type == rhs.type;
  }

  bool operator!=(const Intersection &other)
  {
    return !(*this == other);
  }

  bool isValid()
  {
    return type != NONE;
  }

  string toString()
  {
    stringstream ss;
    ss << "t=" << t << " i=" << i << " b_coords=" << b_coords.toString() << " type=";
    switch (type)
    {
    case SPHERE:
      ss << "SPHERE";
      break;
    case TRIANGLE:
      ss << "TRIANGLE";
      break;
    case LIGHT:
      ss << "LIGHT";
      break;
    case NONE:
      ss << "NONE";
      break;
    }
    return ss.str();
  }
};

const Intersection invalidIntersection(-1, -1, Point::invalidPoint(), NONE);

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b);
void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);

void draw_scene()
{
  for (unsigned int x = 0; x < WIDTH; x++)
  {
    glPointSize(2.0);
    // Do not worry about this usage of OpenGL. This is here just so that we can draw the pixels to the screen,
    // after their R,G,B colors were determined by the ray tracer.
    glBegin(GL_POINTS);
    for (unsigned int y = 0; y < HEIGHT; y++)
    {
      // A simple R,G,B output for testing purposes.
      // Modify these R,G,B colors to the values computed by your ray tracer.
      unsigned char r = (unsigned char)(image[x][y][0] * 255);
      unsigned char g = (unsigned char)(image[x][y][1] * 255);
      unsigned char b = (unsigned char)(image[x][y][2] * 255);

      plot_pixel_display(x, y, r, g, b);
      if (mode == MODE_JPEG)
      {
        if (y < HEIGHT / 2)
          plot_pixel_jpeg(x, y, r, g, b);
        else
          plot_pixel_jpeg(x, y, r, g, b);
      }
    }
    glEnd();
    glFlush();
  }
  printf("Ray tracing completed.\n");
  fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x, y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  // fix: needs to be width, then height here (retarded)
  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in saving\n");
  else
    printf("File saved successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if (strcasecmp(expected, found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parsing error; abnormal program abortion.\n");
    exit(0);
  }
}

void parse_doubles(FILE *file, const char *check, double p[3])
{
  char str[100];
  fscanf(file, "%s", str);
  parse_check(check, str);
  fscanf(file, "%lf %lf %lf", &p[0], &p[1], &p[2]);
  // printf("%s %lf %lf %lf\n", check, p[0], p[1], p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file, "%s", str);
  parse_check("rad:", str);
  fscanf(file, "%lf", r);
  // printf("rad: %f\n", *r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file, "%s", s);
  parse_check("shi:", s);
  fscanf(file, "%lf", shi);
  // printf("shi: %f\n", *shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv, "r");
  if (!file)
  {
    printf("Unable to open input file %s. Program exiting.\n", argv);
    exit(0);
  }

  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file, "%i", &number_of_objects);

  printf("number of objects: %i\n", number_of_objects);

  parse_doubles(file, "amb:", ambient_light);

  for (int i = 0; i < number_of_objects; i++)
  {
    fscanf(file, "%s\n", type);
    printf("%s\n", type);
    if (strcasecmp(type, "triangle") == 0)
    {
      // printf("found triangle\n");
      for (int j = 0; j < 3; j++)
      {
        parse_doubles(file, "pos:", t.v[j].position);
        parse_doubles(file, "nor:", t.v[j].normal);
        parse_doubles(file, "dif:", t.v[j].color_diffuse);
        parse_doubles(file, "spe:", t.v[j].color_specular);
        parse_shi(file, &t.v[j].shininess);
      }

      if (num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if (strcasecmp(type, "sphere") == 0)
    {
      // printf("found sphere\n");

      parse_doubles(file, "pos:", s.position);
      parse_rad(file, &s.radius);
      parse_doubles(file, "dif:", s.color_diffuse);
      parse_doubles(file, "spe:", s.color_specular);
      parse_shi(file, &s.shininess);

      if (num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if (strcasecmp(type, "light") == 0)
    {
      // printf("found light\n");
      parse_doubles(file, "pos:", l.position);
      parse_doubles(file, "col:", l.color);

      if (num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n", type);
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
  glOrtho(0, WIDTH, 0, HEIGHT, 1, -1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0, 0, 0, 0);
  glClear(GL_COLOR_BUFFER_BIT);

  // much faster if done here
  draw_scene();
  if (mode == MODE_JPEG)
    save_jpg();
}

void idle()
{
  // // Hack to make it only draw once.
  // static int once = 0;
  // if (!once)
  // {
  //   draw_scene();
  //   if (mode == MODE_JPEG)
  //     save_jpg();
  // }
  // once = 1;
}

Intersection intersectSphere(Point direction, Point P0, int i, bool shadowRay)
{
  direction.normalize();

  Point C(spheres[i].position);
  double r = spheres[i].radius;

  double b = 2 * dot(direction, P0 - C);
  double c = sq(P0.x - C.x) + sq(P0.y - C.y) + sq(P0.z - C.z) - sq(r);
  double determinant = sq(b) - 4 * c;

  // no real solutions
  if (determinant < 0)
    return invalidIntersection;

  double t0, t1;
  if (determinant == 0)
  {
    t0 = t1 = -b / 2;
  }
  else
  {
    double sqrtDet = sqrt(determinant);
    t0 = (-b - sqrtDet) / 2;
    t1 = (-b + sqrtDet) / 2;
  }

  double t;
  if (t0 > 0 && t1 > 0)
  {
    t = min(t0, t1);
  }
  else if (t0 > 0)
  {
    t = t0;
  }
  else if (t1 > 0)
  {
    t = t1;
  }
  else
  {
    return invalidIntersection;
  }

  return Intersection(t, i, Point(), SPHERE);
}

void projectTriangleTo2D(Point &A, Point &B, Point &C, Point &intersection)
{

  // Compute triangle normal
  Point normal = crossProduct(B - A, C - A);

  // Determine dominant axis of the normal
  if (abs(normal.x) >= abs(normal.y) && abs(normal.x) >= abs(normal.z))
  {
    // Project to YZ plane (discard x)
    std::swap(A.x, A.z);
    std::swap(B.x, B.z);
    std::swap(C.x, C.z);
    std::swap(intersection.x, intersection.z);
  }
  else if (abs(normal.y) >= abs(normal.x) && abs(normal.y) >= abs(normal.z))
  {
    // Project to XZ plane (discard y)
    std::swap(A.y, A.z);
    std::swap(B.y, B.z);
    std::swap(C.y, C.z);
    std::swap(intersection.y, intersection.z);
  }
  // Default: Project to XY plane (discard z)

  // Check for degeneracy after projection
  if ((A.x == B.x && A.y == B.y) || (B.x == C.x && B.y == C.y) || (C.x == A.x && C.y == A.y))
  {
    throw std::runtime_error("Degenerate triangle after projection");
  }

  // // check if triangle is in YZ plane
  // if (equal(A.x, B.x) && equal(B.x, C.x) && equal(A.x, C.x))
  // {
  //   std::swap(A.x, A.z);
  //   std::swap(B.x, B.z);
  //   std::swap(C.x, C.z);
  //   std::swap(intersection.x, intersection.z);
  // }

  // // check if triangle is in XZ plane
  // else if (equal(A.y, B.y) && equal(B.y, C.y) && equal(A.y, C.y))
  // {
  //   std::swap(A.y, A.z);
  //   std::swap(B.y, B.z);
  //   std::swap(C.y, C.z);
  //   std::swap(intersection.y, intersection.z);
  // }

  // // project to 2D
  // A.z = 0;
  // B.z = 0;
  // C.z = 0;
  // intersection.z = 0;
}

Point computeBarycentricCoords(const Point &A, const Point &B, const Point &C, const Point &P)
{
  double area_ABC = (A.x * (B.y - C.y) + B.x * (C.y - A.y) + C.x * (A.y - B.y));
  double area_PBC = (P.x * (B.y - C.y) + B.x * (C.y - P.y) + C.x * (P.y - B.y));
  double area_PCA = (A.x * (P.y - C.y) + P.x * (C.y - A.y) + C.x * (A.y - P.y));
  double area_PAB = (A.x * (B.y - P.y) + B.x * (P.y - A.y) + P.x * (A.y - B.y));

  double alpha = area_PBC / area_ABC;
  double beta = area_PCA / area_ABC;
  double gamma = area_PAB / area_ABC;

  return Point(alpha, beta, gamma);
}

bool checkForIntersection(const Point &c)
{
  return isNormalized(c.x) && isNormalized(c.y) && isNormalized(c.z) && equal(c.x + c.y + c.z, 1.0);
}

Intersection intersectTriangle(Point direction, Point P0, int i, bool shadowRay)
{
  Point A(triangles[i].v[0]);
  Point B(triangles[i].v[1]);
  Point C(triangles[i].v[2]);
  Point normal = crossProduct(B - A, C - A).normalize();

  double d = -dot(normal, A);
  double numerator = -(dot(normal, P0) + d);
  double denominator = dot(normal, direction);

  if (equal(denominator, 0))
    return invalidIntersection;

  double t = numerator / denominator;

  // t <= 0
  if (t < 1e-10)
    return invalidIntersection;

  Point intersectionPoint = P0 + t * direction;

  projectTriangleTo2D(A, B, C, intersectionPoint);
  Point bary_coords = computeBarycentricCoords(A, B, C, intersectionPoint);

  if (checkForIntersection(bary_coords))
    return {t, i, bary_coords, TRIANGLE};
  else
    return invalidIntersection;
}

Point computePhongTriangleNormal(Point b_coords, int triangle_idx)
{
  Triangle triangle = triangles[triangle_idx];
  Point N1(triangle.v[0].normal);
  Point N2(triangle.v[1].normal);
  Point N3(triangle.v[2].normal);
  return (b_coords.x * N1 + b_coords.y * N2 + b_coords.z * N3).normalize();
}

Point computePhongKd(Point b_coords, int triangle_idx)
{
  Triangle triangle = triangles[triangle_idx];
  Point kd1(triangle.v[0].color_diffuse);
  Point kd2(triangle.v[1].color_diffuse);
  Point kd3(triangle.v[2].color_diffuse);
  return (b_coords.x * kd1 + b_coords.y * kd2 + b_coords.z * kd3);
}

Point computePhongKs(Point b_coords, int triangle_idx)
{
  Triangle triangle = triangles[triangle_idx];
  Point ks1(triangle.v[0].color_specular);
  Point ks2(triangle.v[1].color_specular);
  Point ks3(triangle.v[2].color_specular);
  return (b_coords.x * ks1 + b_coords.y * ks2 + b_coords.z * ks3);
}

double computePhongSh(Point b_coords, int triangle_idx)
{
  Triangle triangle = triangles[triangle_idx];
  double sh1 = triangle.v[0].shininess;
  double sh2 = triangle.v[1].shininess;
  double sh3 = triangle.v[2].shininess;
  return dot(b_coords, Point(sh1, sh2, sh3));
}

Point computePhongIllumination(Point lightColor, Point kd, Point L, Point N, Point ks, Point R, Point V, double sh)
{
  double LdotN = dot(L, N);
  double RdotV = dot(R, V);
  if (LdotN < 0)
    LdotN = 0;
  if (RdotV < 0)
    RdotV = 0;

  double r = lightColor.x * (kd.x * LdotN + ks.x * pow(RdotV, sh));
  double g = lightColor.y * (kd.y * LdotN + ks.y * pow(RdotV, sh));
  double b = lightColor.z * (kd.z * LdotN + ks.z * pow(RdotV, sh));

  return Point(r, g, b);
}

Intersection findClosestIntersection(Point ray, Point P0, bool shadowRay)
{
  ray.normalize();

  Intersection closest(std::numeric_limits<double>::max(),
                       -1,
                       Point::invalidPoint(),
                       NONE);

  // calculate intersections with triangles
  for (int i = 0; i < num_triangles; i++)
  {
    Intersection curr = intersectTriangle(ray, P0, i, shadowRay);
    if (curr.isValid() && curr.t < closest.t)
    {
      closest = curr;
    }
  }

  // compute intersections with spheres
  for (int i = 0; i < num_spheres; i++)
  {
    Intersection curr = intersectSphere(ray, P0, i, shadowRay);
    if (curr.isValid() && curr.t < closest.t)
    {
      closest = curr;
    }
  }

  return closest;
}

// returns true if ray is unobstructed, false if blocked
bool shootShadowRay(Point ray, Point P0, int i)
{
  Point light(lights[i].position);
  double t_light = (light.x - P0.x) / ray.x;
  Intersection closest = findClosestIntersection(ray, P0, true);

  // if no intersection, the shadow ray is not blocked
  if (!closest.isValid())
    return true;

  // if there is an intersection, check if it's greater than t_light
  //  if it's greater than t_light, the intersection occurs above/past
  //  the light, so the ray is unobstructed to the light and we return true.
  return closest.t > t_light;
}

/*
 * Shoots out a ray originating at P0. Returns final color
 */
Point shootRay(Point ray, Point P0)
{
  Intersection closest = findClosestIntersection(ray, P0, false);

  Point color(ambient_light);

  // no intersection
  if (closest.type == NONE)
  {
    color = background;
  }
  else
  {
    Point intersection = P0 + closest.t * ray;

    // shoot shadow ray to each light source and add up contributions
    for (int i = 0; i < num_lights; i++)
    {
      Light light = lights[i];
      Point shadowRay = (Point(light.position) - intersection).normalize();

      // fix: make intersection point slightly off the surface to prevent
      // surface from intersecting itself
      Point raisedIntersection = intersection + 2e-5 * shadowRay;

      // shadow ray not blocked
      if (shootShadowRay(shadowRay, raisedIntersection, i))
      {
        // compute phong lighting for each color channel
        Point N(0, 0, 0);
        Point kd(0, 0, 0);
        Point ks(0, 0, 0);
        double sh = 0;
        if (closest.type == TRIANGLE)
        {
          N = computePhongTriangleNormal(closest.b_coords, closest.i);
          kd = computePhongKd(closest.b_coords, closest.i);
          ks = computePhongKs(closest.b_coords, closest.i);
          sh = computePhongSh(closest.b_coords, closest.i);
        }
        else
        {
          N = (intersection - Point(spheres[closest.i].position)).normalize();
          kd = Point(spheres[closest.i].color_diffuse);
          ks = Point(spheres[closest.i].color_specular);
          sh = spheres[closest.i].shininess;
        }

        Point L = shadowRay;
        Point V = (camera - intersection).normalize();
        Point R = 2 * dot(L, N) * N - L;

        // compute r contribution
        color = color + computePhongIllumination(Point(light.color), kd, L, N, ks, R, V, sh);
      }
    }
  }

  color.clampValues();
  return color;
}

void raytrace()
{
  double a = (double)WIDTH / (double)HEIGHT;
  double t = tan(fov * M_PI / 360.0);

  Point topLeft(-a * t, t, -1);

  double width = 2 * a * t;
  double height = 2 * t;
  double x_step = width / (double)WIDTH;
  double y_step = height / (double)HEIGHT;

  cout << "num_triangles: " << num_triangles << endl;
  cout << "num_spheres:   " << num_spheres << endl;
  cout << "num_lights:    " << num_lights << endl;
  cout << endl;

  /**
   * DEBUGGING
   */

  for (unsigned x = 0; x < WIDTH; x++)
  {
    for (unsigned y = 0; y < HEIGHT; y++)
    {
      image[x][HEIGHT - y - 1][0] = 1;
      image[x][HEIGHT - y - 1][1] = 1;
      image[x][HEIGHT - y - 1][2] = 1;
    }
  }

  unsigned x_start = 0;
  unsigned x_end = WIDTH;
  unsigned y_start = 0;
  unsigned y_end = HEIGHT;

  for (unsigned x = x_start; x < x_end; x++)
  {
    for (unsigned y = y_start; y < y_end; y++)
    {
      // cout << "x: " << x << " y: " << y << endl;
      Point ray = topLeft;
      ray.x += x * x_step;
      ray.y -= y * y_step;
      ray.normalize();

      Point finalColor = shootRay(ray, camera);
      image[x][HEIGHT - y - 1][0] = finalColor.x;
      image[x][HEIGHT - y - 1][1] = finalColor.y;
      image[x][HEIGHT - y - 1][2] = finalColor.z;
    }
  }

  int num_pixels = (x_end - x_start) * (y_end - y_start);
  cout << "Drew " << num_pixels << " pixels." << endl;
}

int main(int argc, char **argv)
{
  if ((argc < 2) || (argc > 3))
  {
    printf("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if (argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if (argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc, argv);
  loadScene(argv[1]);

  try
  {
    raytrace();
  }
  catch (std::exception &e)
  {
    cout << e.what() << endl;
  }

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0, 0);
  glutInitWindowSize(WIDTH, HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
#ifdef __APPLE__
  // This is needed on recent Mac OS X versions to correctly display the window.
  glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
#endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
