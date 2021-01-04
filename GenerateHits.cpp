#include "GenerateHits.h"
#define GENERATION
#define Add_Noise
#define x_z_only

void GenerateHits()
{
	
	double xval[nstrips][nplanes];//[nplanes];
	double yval[nstrips][nplanes];//[nplanes];
	double zval[nplanes];
	
	int xhit[nstrips][nplanes];
	int yhit[nstrips][nplanes];
	
	std::default_random_engine generator;
	
	for(int ip=0; ip<nplanes; ip++){
		
		zval[ip] = plane_gap*(ip);
		
		for(int js=0; js<nstrips; js++){
			
			xval[js][ip] = (js)*strip_gap + 0.5*(1+2*js)*strip_width - (0.5*nstrips*strip_width) - (0.5*nstrips-1)*strip_gap - 0.5*strip_gap;
			yval[js][ip] = (js)*strip_gap + 0.5*(1+2*js)*strip_width - (0.5*nstrips*strip_width) - (0.5*nstrips-1)*strip_gap - 0.5*strip_gap;
			
			}
		}
	
	srand (time(NULL));
	
	ofstream fp;
	fp.open("RPC_hits_MC.txt");
	
	for(int iter=0; iter<1000; iter++){
	
//	fp<<"Entry "<<iter+1<<endl;
	
	int xhit_id[nstrips]; 
	int yhit_id[nstrips]; 
	
	
	bool xhit_strips[nstrips][nplanes]={{0}};
	bool yhit_strips[nstrips][nplanes]={{0}};
	
	// generation starts

	#ifdef GENERATION
	
	// add real st. line //
	
	#ifdef x_z_only
	
	// generate c first, then m //

	double cx_true = (double) (rand()%nplanes)*plane_gap ;//((double) rand() / (RAND_MAX))*28;
	double cy_true = (double) (rand()%nplanes)*plane_gap ;//((double) rand() / (RAND_MAX))*28;

//	double theta2 = atan(10.5/max(1.e-6,c_true));
	double theta2 = min(atan(10.5/max(1.e-6,cx_true)),atan(10.5/max(1.e-6,zval[nplanes-1]-cx_true)));
	double theta1 = 0.;//atan(10.5/max(1.e-6,28.-c_true));
	double itrandom = ((double) rand() / (RAND_MAX));
	double theta = acos(1-itrandom*(cos(theta1)-cos(theta2)));
	double phi =  M_PI * double(rand()%2);
	theta *= cos(phi);
	double mx_true = 1./(tan(theta));//*cos(phi));
	double my_true = 1./(tan(theta)*sin(phi));
	
	// end //
	
	#else
	
	// generate theta with x,z constraints //
	
	double x_true = (((double) rand() / (RAND_MAX)) * (xval[nstrips-1][nplanes-1] + strip_width - xval[0][nplanes-1])) + (xval[0][nplanes-1] - 0.5*strip_width);
	double y_true = (((double) rand() / (RAND_MAX)) * (yval[nstrips-1][nplanes-1] + strip_width - yval[0][nplanes-1])) + (yval[0][nplanes-1] - 0.5*strip_width);

	double theta2 = atan((xval[nstrips-1][nplanes-1] - x_true)*1./(zval[nplanes-1]-zval[0]));
	double theta1 = atan((x_true - xval[0][nplanes-1])*1./(zval[nplanes-1]-zval[0]));
	double itrandom = ((double) rand() / (RAND_MAX));
	double theta = acos(1-itrandom*(cos(theta1)-cos(theta2)));
	double phi = 2*M_PI*((double) rand() / (RAND_MAX));
	double mx_true = 1./(tan(theta)*cos(phi));
	double my_true = 1./(tan(theta)*sin(phi));
	double cx_true = (zval[nplanes-1] - mx_true*x_true) ;
	double cy_true = (zval[nplanes-1] - my_true*y_true) ;
	//end
	
	#endif
	/*
	cout<<"true theta "<<theta<<"("<<theta*180*1./3.14<<")"<<endl;
	cout<<"true phi "<<phi<<endl;
	cout<<"true m "<<mx_true<<endl;
	cout<<"cx_true "<<cx_true<<" cy_true "<<cy_true<<endl;
	cout<<"cxs_true "<<cx_true*sin(theta)<<endl;
	*/
	int nxtrue[nplanes] = {0};
	int nytrue[nplanes] = {0};
	
	for(int ip=0; ip<nplanes; ip++){
		
		double tryx = (zval[ip] - cx_true)*1./mx_true;
		
		for(int ix=0; ix<nstrips; ix++){
			if(tryx >= (xval[ix][ip] - 0.5*strip_width) && tryx <= (xval[ix][ip] + 0.5*strip_width)){
				xhit_strips[ix][ip] = 1;
				double ispread = ((double) rand() / (RAND_MAX));//rand()%10;
				if(ispread>= 0 && ispread<0.05) { 
					xhit_strips[ix][ip] = 0; 
					}
				if(ispread>= 0.05 && ispread<0.15) {
					if(ix>0 && ix<(nstrips-1)){
						xhit_strips[ix-1][ip] = xhit_strips[ix+1][ip] = 1;
						}
					if(ix==0){
						xhit_strips[ix+1][ip] = 1;
						}
					if(ix==(nstrips-1)){
						xhit_strips[ix-1][ip] = 1;
						}	
				}
				if(ispread>= 0.75 && ispread<1.) {
					if(ix>0 && ix<(nstrips-1)){
						if(tryx <= xval[ix][ip]) { xhit_strips[ix-1][ip] = 1; }
						if(tryx > xval[ix][ip]) { xhit_strips[ix+1][ip] = 1; }
					}
					if(ix==0){
						if(tryx > xval[ix][ip]) { xhit_strips[ix+1][ip] = 1; }
						}
					if(ix==(nstrips-1)){
						if(tryx <= xval[ix][ip]) { xhit_strips[ix-1][ip] = 1; }
					}	
				}
				nxtrue[ip]++;
				}
		}
		
		double tryy = (zval[ip] - cx_true)*1./my_true;
		for(int iy=0; iy<nstrips; iy++){
			if(tryy >= (yval[iy][ip] - 0.5*strip_width) && tryy <= (yval[iy][ip] + 0.5*strip_width)){
				yhit_strips[iy][ip] = 1;
				double ispread = ((double) rand() / (RAND_MAX));//rand()%10;
				if(ispread>= 0 && ispread<0.05) { yhit_strips[iy][ip] = 0; }
				if(ispread>= 0.05 && ispread<0.15) {
					if(iy>0 && iy<(nstrips-1)){
						yhit_strips[iy-1][ip] = yhit_strips[iy+1][ip] = 1;
						}
					if(iy==0){
						yhit_strips[iy+1][ip] = 1;
						}
					if(iy==(nstrips-1)){
						yhit_strips[iy-1][ip] = 1;
						}	
				}
				if(ispread>= 0.75 && ispread<1.) {
					if(iy>0 && iy<(nstrips-1)){
						if(tryy <= yval[iy][ip]) { yhit_strips[iy-1][ip] = 1; }
						if(tryy > yval[iy][ip]) { yhit_strips[iy+1][ip] = 1; }
					}
					if(iy==0){
						if(tryy > yval[iy][ip]) { yhit_strips[iy+1][ip] = 1; }
						}
					if(iy==(nstrips-1)){
						if(tryy <= yval[iy][ip]) { yhit_strips[iy-1][ip] = 1; }
					}	
				}
				nytrue[ip]++;
				}
		}	
	}
	
	// real st. line done //
	
//	cout<<"RAW hits (signal):"<<endl;
	#ifndef Add_Noise
	
	for(int ip=0; ip<nplanes; ip++){
		for(int ix1=0; ix1<nstrips; ix1++){
			fp<<xhit_strips[ix1][ip]<<"\t";
		}
		fp<<endl;
	}
	fp<<endl
	
	#endif
	
//	if(!(ntrue[0]>0 && ntrue[nplanes-1]>0)) return;
//	if(ntrue_total<4) return;

	double proj_x2 = (zval[nplanes-1] - cx_true)*1./mx_true;
	double proj_x1 = (zval[0] - cx_true)*1./mx_true;
	if(proj_x2 > xval[nstrips-1][nplanes-1] || proj_x1 < xval[0][0]) return;
	if(proj_x1 > xval[nstrips-1][nplanes-1] || proj_x2 < xval[0][0]) return;
	
	double proj_y2 = (zval[nplanes-1] - cy_true)*1./my_true;
	double proj_y1 = (zval[0] - cy_true)*1./my_true;
	if(proj_y2 > yval[nstrips-1][nplanes-1] || proj_y1 < yval[0][0]) return;
	if(proj_y1 > yval[nstrips-1][nplanes-1] || proj_y2 < yval[0][0]) return;
	
	#ifdef Add_Noise
	
	int max_noise_x = 2;
	int max_noise_y = 2;
	
	for(int ip=0; ip<nplanes; ip++){
		
		if((rand()%2)==0) continue;
		
		int nxhits = rand()%(max_noise_x+1);//(nstrips+1);
		gethits(xhit_id,nxhits);
		
		for(int ixs=0; ixs<nxhits; ixs++){
			xhit_strips[xhit_id[ixs]-1][ip] = 1;		
			}
		
		int nyhits = rand()%(max_noise_y+1);//(nstrips+1);
		gethits(yhit_id,nyhits);
		
		for(int iys=0; iys<nyhits; iys++){
			yhit_strips[yhit_id[iys]-1][ip] = 1;
			}
	
	}
	
	#endif
	// noise done //
	
//	cout<<"RAW hits (signal+noise):"<<endl;
	
	for(int ip=0; ip<nplanes; ip++){
		for(int ix1=0; ix1<nstrips; ix1++){
			fp<<xhit_strips[ix1][ip]<<"\t";
		}
		fp<<endl;
	}
//	fp<<endl;
	
	#endif
	
	// generation ends
	}//iter

	fp.close();
	
}

