#include <stdio.h> 
#include <stdlib.h> 
#include <time.h>
#include <iostream>
#include <math.h>
#include <vector>
//#include <random>

static const int nstrips = 8;
static const int nplanes = 8;

const double strip_gap = 0.2;
const double strip_width = 2.8;
const double plane_gap = 4;

struct fitparam{
  float x1;
  float x2;
  float z1;
  float z2;
  float theta;
  float rho;
};

struct point{
 float x;
 float z;
};

void trial();
