#include <stdio.h> 
#include <stdlib.h> 
#include <time.h>
#include <iostream>
#include <fstream>
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

void gethits(int *dst, int nhits)
{

int sz, pos, src[nstrips];

for (int i = 0; i < sizeof(src)/sizeof(*src); i++){
	src[i] = i + 1;
}

sz = nstrips;
for (int i = 0; i < nhits; i++) {
	pos = rand() % sz;
	dst[i] = src[pos];
	src[pos] = src[sz-1];
	sz--;
}
	
};

/*
double chi2_func(double yval, double xval, double m, double c){
	
double ypred = (m*xval+c);
double yerr = m*strip_width*1./sqrt(12);
double chi2 = (yval - ypred)*(yval - ypred)*1./(yerr*yerr);
return chi2;
};
*/
void Cosmic_trig();
