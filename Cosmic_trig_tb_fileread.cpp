#include "Cosmic_trig.h"

#define Remove_Noise
#define Optimize
#define Optimize_Loop
#define x_z_only

#define Break_chain

using namespace std;

int main(void){

#ifndef __SYNTHESIS__
ifstream fp("/home/suman/INO/INO_Data.txt");
if(!(fp.is_open())) return 0;

string line;
#endif

float zval[nplanes];

for(int ip=0; ip<nplanes; ip++){
	zval[ip] = plane_gap*(ip);
}
	
for(int iter=0; iter<4; iter++){
	
//	cout<<"======== Iter =========== "<<iter+1<<endl;
	
	#ifndef __SYNTHESIS__
	if(fp.eof()){  break; }
	
	if(iter==0){
	getline(fp,line);	
	cout<<line<<endl;
	}
	if(iter>0){
	fp.ignore();
	getline(fp,line);	
	cout<<line<<endl;
	}
	#endif

	bool xhitx_strips[nstrips][nplanes]={{0}};
	bool yhitx_strips[nstrips][nplanes]={{0}};
	
	/*
	bool xhitx_strips[nstrips][nplanes]=
	{{0,1,1,1,0,0,1,0},	
	 {0,0,1,1,0,0,0,0},	
	 {0,0,1,1,0,0,1,0},	
	 {0,0,1,1,1,0,0,0},	
	 {1,0,0,1,0,0,0,0},	
	 {0,0,0,1,1,0,0,0},	
	 {0,0,1,1,1,0,0,0},	
	 {0,0,0,0,0,0,0,0}};
	 
	 bool yhitx_strips[nstrips][nplanes]=
	{{0,0,0,0,0,1,0,0},	
	 {0,0,0,0,0,1,0,0},	
	 {0,0,0,0,0,1,0,0},	
	 {0,0,0,0,0,1,0,0},	
	 {1,0,0,1,1,0,0,0},	
	 {0,0,0,0,1,0,0,0},	
	 {0,0,0,0,1,0,0,0},	
	 {0,0,0,0,0,0,1,1}};
	 */

	#ifndef __SYNTHESIS__

	for(int ip=0; ip<nplanes; ip++){
		for(int ix1=0; ix1<nstrips; ix1++){
			int xx;
			fp>>xx;
			xhitx_strips[ix1][ip] = (xx>0)?true:false;
		}
	}
	
	for(int ip=0; ip<nplanes; ip++){
		for(int ix1=0; ix1<nstrips; ix1++){
			int xx;
			fp>>xx;
			yhitx_strips[ix1][ip] = (xx>0)?true:false;
		}
	}

	#endif
	
	bool xzpass;
	
	Cosmic_trig_top(xhitx_strips,yhitx_strips,zval,xzpass);
	
	cout<<"Event pass: "<<xzpass<<endl;
	
}//iter

#ifndef __SYNTHESIS__
  fp.close();
#endif

return(0);

}
