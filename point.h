#include <sstream>
#include <string>
#include <iostream>
#include <stdexcept>
#include "objects.h"
#include "utils.h"

#ifndef POINT_H
#define POINT_H

using namespace std;

// Represents one spline control point.
struct Point
{
  double x, y, z;

  Point() : valid(false) {}

  Point(double x, double y, double z) : x(x), y(y), z(z), valid(true) {}

  Point(double *point) : x(point[0]), y(point[1]), z(point[2]), valid(true) {}

  Point(const Vertex &v) : x(v.position[0]), y(v.position[1]), z(v.position[2]), valid(true) {}

  string toString() const
  {
    checkIfValid();
    stringstream ss;
    ss << "(" << x << ", " << y << ", " << z << ")";
    return ss.str();
  }

  Point &print()
  {
    checkIfValid();
    cout << toString() << endl;
    return *this;
  }

  Point operator+(const Point &rhs) const
  {
    checkIfValid();
    return Point(x + rhs.x, y + rhs.y, z + rhs.z);
  }

  Point operator-() const
  {
    checkIfValid();
    return Point(-x, -y, -z);
  }

  Point operator-(const Point &rhs) const
  {
    checkIfValid();
    return *this + -rhs;
  }

  bool operator==(const Point &other) const
  {
    checkIfValid();
    return x == other.x && y == other.y && z == other.z && valid == other.valid;
  }

  Point &normalize()
  {
    checkIfValid();
    double m = magnitude();
    x /= m;
    y /= m;
    z /= m;
    return *this;
  }

  Point& clampValues()
  {
    x = clamp(x);
    y = clamp(y);
    z = clamp(z);
  }

  double magnitude() const
  {
    checkIfValid();
    return sqrt((x * x) + (y * y) + (z * z));
  }

  bool isValid() const
  {
    return valid;
  }

  static Point invalidPoint()
  {
    Point p(0, 0, 0);
    p.valid = false;
    return p;
  }

private:
  bool valid;
  void checkIfValid() const
  {
    if (!valid)
      throw std::invalid_argument("Point is invalid");
  }
};

Point normalize(const Point &p)
{
  return Point(p).normalize();
}

Point crossProduct(const Point &A, const Point &B)
{
  double Cx = A.y * B.z - A.z * B.y;
  double Cy = A.z * B.x - A.x * B.z;
  double Cz = A.x * B.y - A.y * B.x;
  return Point(Cx, Cy, Cz);
}

double dot(const Point &A, const Point &B)
{
  return (A.x * B.x) + (A.y * B.y) + (A.z * B.z);
}

Point operator*(const float a, const Point &p)
{
  return Point(a * p.x, a * p.y, a * p.z);
}

std::ostream &operator<<(std::ostream &os, const Point &p)
{
  os << p.toString();
  return os;
}

#endif