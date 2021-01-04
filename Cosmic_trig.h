#include <stdio.h>
#include <stdlib.h> 
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>
//#include <random>

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

const int nstrips = 8;
const int nplanes = 8;

const float strip_gap = 0.2;
const float strip_width = 2.8;
const float plane_gap = 4;

void Cosmic_trig_top(const bool x[nstrips][nplanes],const bool y[nstrips][nplanes],const float zval[nplanes], bool &xzpass);
