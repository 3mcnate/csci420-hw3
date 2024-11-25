#ifndef UTILS_H
#define UTILS_H

#include <cmath>

float clamp(double x)
{
  if (x < 0)
    return 0;
  if (x > 1)
    return 1;
  return x;
}

double sq(double x)
{
  return x * x;
}

bool equal(double x, double y)
{
  return abs(x - y) < 0.0001;
}

#endif