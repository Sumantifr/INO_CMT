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
	
	bool xhitx_strips[nstrips][nplanes]={{0}};
	
	xhitx_strips[0][0] = 1;
	xhitx_strips[1][0] = 1;
	xhitx_strips[0][1] = 1;
	xhitx_strips[0][2] = 1;
	xhitx_strips[1][2] = 1;
	xhitx_strips[0][3] = 1;
	xhitx_strips[1][3] = 1;
	xhitx_strips[1][4] = 1;
	xhitx_strips[1][5] = 1;

	 
	bool yhitx_strips[nstrips][nplanes]={{0}};
	
	bool xzpass;
	
	Cosmic_trig_top(xhitx_strips,yhitx_strips,zval,xzpass);
	
	cout<<"Event pass: "<<xzpass<<endl;

return(0);

}
