#include "Cosmic_trig.h"

#define Remove_Noise
#define Optimize
#define Optimize_Loop
#define x_z_only

#define Break_chain

using namespace std;

int main(void){

float zval[nplanes];

for(int ip=0; ip<nplanes; ip++){
	zval[ip] = plane_gap*(ip);
}
	
//	cout<<"======== Iter =========== "<<iter+1<<endl;
	
	bool xhitx_strips[nstrips][nplanes]=
	{
	{1,1,0,0,0,0,0,0},
	{1,0,0,0,0,0,0,0},
	{1,1,0,0,0,0,0,0},
	{1,1,0,0,0,0,0,0},
	{0,1,0,0,0,0,0,0},
	{0,1,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0}
	};

	 
	 bool yhitx_strips[nstrips][nplanes]=
	{
	{0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0},
	{0, 0, 0, 0, 0, 0, 0, 0}
	};
	
	bool xzpass;
	
	Cosmic_trig_top(xhitx_strips,yhitx_strips,zval,xzpass);
	
	cout<<"Event pass: "<<xzpass<<endl;

return(0);

}
