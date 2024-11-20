#ifndef UTILS_H
#define UTILS_H

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

#endif