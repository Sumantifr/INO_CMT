#include "trial.h"

using namespace std;
/*
#define GENERATION
#define Add_Noise
#define Remove_Noise
#define Remove_Noise_in_Signal
//#define PRINT_Details
#define Optimize
#define Optimize_Loop
#define x_z_only
*/

void trial()
{
	
	double xval[nstrips][nplanes];//[nplanes];
	double yval[nstrips][nplanes];//[nplanes];
	double zval[nplanes];
	
	double xval_true[nplanes];
	double xval_sig[nplanes][nplanes];
	
	int xhit[nstrips][nplanes];
	int yhit[nstrips][nplanes];
	
//	std::default_random_engine generator;
	
	for(int ip=0; ip<nplanes; ip++){
		
		zval[ip] = plane_gap*(ip);
		
		for(int js=0; js<nstrips; js++){
			
			xval[js][ip] = (js)*strip_gap + 0.5*(1+2*js)*strip_width - (0.5*nstrips*strip_width) - (0.5*nstrips-1)*strip_gap - 0.5*strip_gap;
			yval[js][ip] = (js)*strip_gap + 0.5*(1+2*js)*strip_width - (0.5*nstrips*strip_width) - (0.5*nstrips-1)*strip_gap - 0.5*strip_gap;
			
			xval_sig[js][ip] = xval[js][ip];
			if(ip==0) { cout<<js+1<<" x: "<<xval[js][ip]<<endl; }
			
			}
		
		 cout<<ip+1<<" z: "<<zval[ip]<<endl; 
		
	}
	
//	srand (time(NULL));
	
	int xhit_id[nstrips]; 
	int yhit_id[nstrips]; 
	
	// add noise //
	
	bool xhit_strips[nstrips][nplanes]={{0}};
	bool yhit_strips[nstrips][nplanes]={{0}};
	
	bool xhit_strips_sig[nstrips][nplanes]={{0}};
	
	// generation starts

	#ifdef GENERATION
	
	// add real st. line //
	
	// generate c first, then m //

	double cx_true = 10;//(double) (rand()%nplanes)*plane_gap ;//((double) rand() / (RAND_MAX))*28;
	double cy_true = 10;//(double) (rand()%nplanes)*plane_gap ;//((double) rand() / (RAND_MAX))*28;

//	double theta2 = atan(10.5/max(1.e-6,c_true));
	double theta2 = min(atan(10.5/max(1.e-6,cx_true)),atan(10.5/max(1.e-6,zval[nplanes-1]-cx_true)));
	double theta1 = 0.;//atan(10.5/max(1.e-6,28.-c_true));
	double itrandom = 0.2;//((double) rand() / (RAND_MAX));
	double theta = acos(1-itrandom*(cos(theta1)-cos(theta2)));
	double phi =  M_PI * 0;//double(rand()%2);
	theta *= cos(phi);
	double mx_true = 1./(tan(theta));//*cos(phi));
	double my_true = 1./(tan(theta)*sin(phi));
	
	// end //
	
	
	// cout<<"true theta "<<theta<<"("<<theta*180*1./3.14<<")"<<endl;
	// cout<<"true phi "<<phi<<endl;
	// cout<<"true m "<<mx_true<<endl;
	// cout<<"cx_true "<<cx_true<<" cy_true "<<cy_true<<endl;
	// cout<<"cxs_true "<<cx_true*sin(theta)<<endl;
	
	int nxtrue[nplanes] = {0};
	int nytrue[nplanes] = {0};
	
	for(int ip=0; ip<nplanes; ip++){
		
		double tryx = (zval[ip] - cx_true)*1./mx_true;
		xval_true[ip] = tryx;
		
//		// cout<<"plane:"<<ip+1<<endl;
		for(int ix=0; ix<nstrips; ix++){
//			// cout<<"x,xmin,xmax "<<tryx<<" "<<(xval[ix] - 0.5*strip_width)<<" "<<(xval[ix] + 0.5*strip_width)<<endl;
			if(tryx >= (xval[ix][ip] - 0.5*strip_width) && tryx <= (xval[ix][ip] + 0.5*strip_width)){
				xhit_strips[ix][ip] = 1;
				xhit_strips_sig[ix][ip] = 1;
				double ispread = 0.3;//((double) rand() / (RAND_MAX));//rand()%10;
				if(ispread>= 0 && ispread<0.05) { 
					xhit_strips[ix][ip] = 0; 
					xhit_strips_sig[ix][ip] = 0;
					}
				if(ispread>= 0.05 && ispread<0.15) {
					if(ix>0 && ix<(nstrips-1)){
						xhit_strips[ix-1][ip] = xhit_strips[ix+1][ip] = 1;
						xhit_strips_sig[ix-1][ip] = xhit_strips_sig[ix+1][ip] = 1;
						}
					if(ix==0){
						xhit_strips[ix+1][ip] = 1;
						xhit_strips_sig[ix+1][ip] = 1;
						}
					if(ix==(nstrips-1)){
						xhit_strips[ix-1][ip] = 1;
						xhit_strips_sig[ix-1][ip] = 1;
						}	
				}
				if(ispread>= 0.75 && ispread<1.) {
					if(ix>0 && ix<(nstrips-1)){
						if(tryx <= xval[ix][ip]) { xhit_strips[ix-1][ip] = 1; xhit_strips_sig[ix-1][ip] = 1; }
						if(tryx > xval[ix][ip]) { xhit_strips[ix+1][ip] = 1; xhit_strips_sig[ix+1][ip] = 1; }
					}
					if(ix==0){
						if(tryx > xval[ix][ip]) { xhit_strips[ix+1][ip] = 1; xhit_strips_sig[ix+1][ip] = 1; }
						}
					if(ix==(nstrips-1)){
						if(tryx <= xval[ix][ip]) { xhit_strips[ix-1][ip] = 1; xhit_strips_sig[ix-1][ip] = 1; }
					}	
				}
//				// cout<<ix<<","<<ip<<endl;
				nxtrue[ip]++;
				}
		}
		
		double tryy = (zval[ip] - cx_true)*1./my_true;
//		// cout<<"plane:"<<ip+1<<endl;
		for(int iy=0; iy<nstrips; iy++){
//			// cout<<"x,xmin,xmax "<<tryx<<" "<<(xval[ix] - 0.5*strip_width)<<" "<<(xval[ix] + 0.5*strip_width)<<endl;
			if(tryy >= (yval[iy][ip] - 0.5*strip_width) && tryy <= (yval[iy][ip] + 0.5*strip_width)){
				yhit_strips[iy][ip] = 1;
				double ispread = 0.4;//((double) rand() / (RAND_MAX));//rand()%10;
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
//				// cout<<ix<<","<<ip<<endl;
				nytrue[ip]++;
				}
		}	
	}
	
	// real st. line done //
	
	// cout<<"RAW hits (signal):"<<endl;
	
	for(int ip=0; ip<nplanes; ip++){
		for(int ix1=0; ix1<nstrips; ix1++){
			// cout<<xhit_strips[ix1][ip]<<"\t";
		}
		// cout<<endl;
	}
	
//	if(!(ntrue[0]>0 && ntrue[nplanes-1]>0)) return;
//	if(ntrue_total<4) return;

	double proj_x2 = (zval[nplanes-1] - cx_true)*1./mx_true;
	double proj_x1 = (zval[0] - cx_true)*1./mx_true;
//	// cout<<"proj_x2 "<<proj_x2<<" proj_x1 "<<proj_x1<<endl;
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
		
	//	if((rand()%2)==0) continue;
		
		int nxhits = 2;//rand()%(max_noise_x+1);//(nstrips+1);
	//	// cout<<"nxhits "<<nxhits<<endl;
	//	gethits(xhit_id,nxhits);
		
		// get hits in x-z plane //
		
		int sxz, posx, srcx[nstrips];

		for (int i = 0; i < sizeof(srcx)/sizeof(*srcx); i++){
			srcx[i] = i + 1;
		}

		sxz = nstrips;
		for (int i = 0; i < nxhits; i++) {
			posx = 5 ;//rand() % sxz;
			xhit_id[i] = srcx[posx];
			srcx[posx] = srcx[sxz-1];
			sxz--;
		}
		
		// hits are obtained
		
		for(int ixs=0; ixs<nxhits; ixs++){
			xhit_strips[xhit_id[ixs]-1][ip] = 1;		
			}
		
		int nyhits = 2;//rand()%(max_noise_y+1);//(nstrips+1);
	//	// cout<<"nyhits "<<nyhits<<endl;
	//	gethits(yhit_id,nyhits);
		
		// get hits in y-z plane
		
		int syz, posy, srcy[nstrips];

		for (int i = 0; i < sizeof(srcy)/sizeof(*srcy); i++){
			srcy[i] = i + 1;
		}

		syz = nstrips;
		for (int i = 0; i < nyhits; i++) {
			posy = 5;//rand() % syz;
			yhit_id[i] = srcy[posy];
			srcy[posy] = srcy[syz-1];
			syz--;
		}
		
		// hits are obtained
		
		for(int iys=0; iys<nyhits; iys++){
			yhit_strips[yhit_id[iys]-1][ip] = 1;
			}
	
	}
	
	#endif
	// noise done //
	
//	 cout<<"RAW hits (signal+noise):"<<endl;
	
	for(int ip=0; ip<nplanes; ip++){
		for(int ix1=0; ix1<nstrips; ix1++){
	//		 cout<<xhit_strips[ix1][ip]<<"\t";
		}
	//	 cout<<endl;
	}
	
	#endif
	
	// generation ends
	
	
	#ifdef Remove_Noise
	
	for(int ip=0; ip<nplanes; ip++){
		
		// if there are 4 strips with hits in a single layer, discard the layer //
		
		int nstr_hit = 0; 
		for(int ix1=0; ix1<nstrips; ix1++){
			if(xhit_strips[ix1][ip]){
				nstr_hit++;
			}
		}
		if(nstr_hit>=5){
			for(int ix1=0; ix1<nstrips; ix1++){
				xhit_strips[ix1][ip] = 0;
				}
			}
		
		nstr_hit = 0;
		for(int iy1=0; iy1<nstrips; iy1++){
			if(yhit_strips[iy1][ip]){
				nstr_hit++;
			}
		}
		if(nstr_hit>=5){
			for(int iy1=0; iy1<nstrips; iy1++){
				yhit_strips[iy1][ip] = 0;
				}
			}
		
			
		// done//
		
		// if there are hits in 3 consecutive strips, merge them //
		
		for(int ix1=0; ix1<(nstrips-2); ix1++){
			if(xhit_strips[ix1][ip] && xhit_strips[ix1+1][ip] && xhit_strips[ix1+2][ip]){
				xhit_strips[ix1][ip] = xhit_strips[ix1+2][ip] = 0;
				}
		}
		
		for(int iy1=0; iy1<(nstrips-2); iy1++){
			if(yhit_strips[iy1][ip] && yhit_strips[iy1+1][ip] && yhit_strips[iy1+2][ip]){
				yhit_strips[iy1][ip] = yhit_strips[iy1+2][ip] = 0;
				}
		}
		
		//  if there are hits in 2 consecutive strips, merge them //
		
		for(int ix1=0; ix1<(nstrips-1); ix1++){
			if(xhit_strips[ix1][ip] && xhit_strips[ix1+1][ip]){
				xhit_strips[ix1+1][ip] = 0;
				xval[ix1][ip] = 0.5*(xval[ix1][ip]+xval[ix1+1][ip]);
//				// cout<<"ix1 "<<ix1+1<<" xval "<<xval[ix1][ip]<<endl;
				}
		}
		
		for(int iy1=0; iy1<(nstrips-1); iy1++){
			if(yhit_strips[iy1][ip] && yhit_strips[iy1+1][ip]){
				yhit_strips[iy1+1][ip] = 0;
				yval[iy1][ip] = 0.5*(yval[iy1][ip]+yval[iy1+1][ip]);
				}
		}
		
		// done
	}//
	
	#endif
	
	// cout<<"Skimmed hits:"<<endl;
	
	for(int ip=0; ip<nplanes; ip++){
		for(int ix1=0; ix1<nstrips; ix1++){
			// cout<<xhit_strips[ix1][ip]<<"\t";
		}
		// cout<<endl;
	}
	
	
	const int theta_bins = 10;//40;
//	double theta_max = M_PI + 0.05; double theta_min = 0 - 0.05;
	double theta_max = 0.5*M_PI + 0.05; double theta_min = -0.5*M_PI;
	double theta_binsize = (theta_max-theta_min)*1./theta_bins;
	double theta1_max = 0.5*M_PI + 0.05; double theta1_min = -0.5*M_PI - 0.05;
	double theta1_binsize = (theta1_max-theta1_min)*1./theta_bins;
	
	const int rho_bins = 15;//30;
	double rho_max = 28.+ 1; double rho_min = -28-1;
	double rho_binsize = (rho_max-rho_min)*1./rho_bins;
	
	const int cs_bins = 60;//(int(zval[nplanes-1]-zval[0]));
	double cs_max = 28. + 1; double cs_min = -28-1;
	double cs_binsize = (cs_max-cs_min)*1./cs_bins;
	
	const int ntfbins = 500;
	double tfbinsize = (theta1_max-theta1_min)*1./ntfbins;
	
	const int ncfbins = 300;//(int(zval[nplanes-1]-zval[0]));
	double cfbinsize = (cs_max-cs_min)*1./ncfbins;
	
	static const int max_entry = 100;
	
	
	int boxes[theta_bins][rho_bins] = {{0}};
	int boyes[theta_bins][rho_bins] = {{0}};
	int toxes[theta_bins][cs_bins] = {{0}};
	
	
	// fill grids in x-z plane
	
	for(int ip=0; ip<nplanes; ip++){
		
		for(int ix1=0; ix1<nstrips; ix1++){
			
			if(xhit_strips[ix1][ip]){
				
				for(int jp=(ip+1); jp<nplanes; jp++){
				
//					if(jp>(ip+3)) continue; // forming pairs only between 3 consecutive layers
				
					for(int ix2=0; ix2<nstrips; ix2++){
					
						if(xhit_strips[ix2][jp]){
							
							double theta_test, tantheta_test, rho_test;
							
							tantheta_test = -1.*(xval[ix2][jp] - xval[ix1][ip])*1./(zval[jp] - zval[ip]);
							theta_test = atan(tantheta_test);
							rho_test = (xval[ix1][ip] + tantheta_test * zval[ip])*1./sqrt(1+tantheta_test*tantheta_test);
							rho_test += (xval[ix2][jp] + tantheta_test * zval[jp])*1./sqrt(1+tantheta_test*tantheta_test);
							rho_test = 0.5* rho_test;
//							// cout<<"(ix1,iz1):"<<" ("<<ix1+1<<","<<ip+1<<")"<<"(ix2,iz2):"<<" ("<<ix2+1<<","<<jp+1<<")"<<" (theta,rho): ("<<theta_test<<","<<rho_test<<")\n";
							
							if(rho_test>=rho_min && rho_test<rho_max && theta_test>=theta_min && theta_test<theta_max){
								
								int itheta_bin = int((theta_test-theta_min)*1/theta_binsize);
								int irho_bin =	int((rho_test-rho_min)*1./rho_binsize);
								boxes[itheta_bin][irho_bin] += 1;
								
							}
							
							
						}
					
					}
				}
				
			}
			
		}
	}
	
	// fill grids in y-z plane
	
	for(int ip=0; ip<nplanes; ip++){
		
		for(int iy1=0; iy1<nstrips; iy1++){
			
			if(yhit_strips[iy1][ip]){
				
				for(int jp=(ip+1); jp<nplanes; jp++){
				
//					if(jp>(ip+3)) continue; // forming pairs only between 3 consecutive layers
				
					for(int iy2=0; iy2<nstrips; iy2++){
					
						if(yhit_strips[iy2][jp]){
							
							double theta_test, tantheta_test, rho_test;
							
							tantheta_test = -1.*(yval[iy2][jp] - yval[iy1][ip])*1./(zval[jp] - zval[ip]);
							theta_test = atan(tantheta_test);
							rho_test = (yval[iy1][ip] + tantheta_test * zval[ip])*1./sqrt(1+tantheta_test*tantheta_test);
							rho_test += (yval[iy2][jp] + tantheta_test * zval[jp])*1./sqrt(1+tantheta_test*tantheta_test);
							rho_test = 0.5* rho_test;
//							// cout<<"(ix1,iz1):"<<" ("<<ix1+1<<","<<ip+1<<")"<<"(ix2,iz2):"<<" ("<<ix2+1<<","<<jp+1<<")"<<" (theta,rho): ("<<theta_test<<","<<rho_test<<")\n";
							
							if(rho_test>=rho_min && rho_test<rho_max && theta_test>=theta_min && theta_test<theta_max){
								
								int itheta_bin = int((theta_test-theta_min)*1/theta_binsize);
								int irho_bin =	int((rho_test-rho_min)*1./rho_binsize);
								boyes[itheta_bin][irho_bin] += 1;
								
							}
							
						}
					
					}
				}
				
			}
			
		}
	}

	int xxx=0;
	
	vector<int> final_xpars;
	vector<int> final_ypars;

	int it_max = -1; int ir_max = -1; 
	int ity_max = -1; int iry_max = -1;
	
	double max_con = -100;
	double max_ycon = -100;
	
	for(int it=0; it<theta_bins; it++){
		for(int jc=0; jc<rho_bins; jc++){
			if(boxes[it][jc]>= max_con){
				max_con = boxes[it][jc];
				it_max = it;
				ir_max = jc;
				}
			}
		}
		
	for(int it=0; it<theta_bins; it++){
		for(int jc=0; jc<rho_bins; jc++){
			if(boyes[it][jc]>= max_ycon){
				max_ycon = boyes[it][jc];
				ity_max = it;
				iry_max = jc;
				}
			}
		}	
	
	for(int ic=0; ic<boxes[it_max][ir_max]; ic++){
		
	//	final_xpars.push_back(1);
	
		}
		
	for(int ic=0; ic<boyes[ity_max][iry_max]; ic++){
		
	//	final_ypars.push_back(1);
	
		}	
	
	double t_opt = 	(theta_min+theta_binsize*(it_max+0.5));
	double r_opt = 	(rho_min+rho_binsize*(ir_max+0.5));

	double ty_opt = 	(theta_min+theta_binsize*(ity_max+0.5));
	double ry_opt = 	(rho_min+rho_binsize*(iry_max+0.5));
	/*
	double deno, t_s, r_s;
	double deno1, ty_s, ry_s;
	
	t_s = -100; r_s =  -100;
	ty_s = -100; ry_s =  -100;
//	t_s = t_opt; r_s = r_opt;
	
	#ifdef Optimize_Loop
	 
	const int nloop = 5;
	
 // loop in x-z  plane
	
	for(int iloop=0; iloop<nloop; iloop++){
	 
		if(fabs(t_s - t_opt) > 2*theta_binsize || fabs(r_s - r_opt) > 2*rho_binsize){
		 
		 if(iloop > 0){
		 
			for(int it_1=0; it_1<theta_bins; it_1++){
				for(int jc_1=0; jc_1<rho_bins; jc_1++){
		 
						if( (t_s > (theta_min + (it_1)*theta_binsize)) && (t_s < (theta_min + (it_1+1)*theta_binsize))
							&& (r_s > (rho_min + (jc_1)*rho_binsize)) && (r_s < (rho_min + (jc_1+1)*rho_binsize)) )
						{	
							it_max = it_1; ir_max = jc_1;
						}
					}
				}
			
				t_opt = t_s;
				r_opt = r_s;
			
			}
			
		 deno = 0;
		 t_s =0; r_s = 0;
		 
		 if(it_max>=0 && it_max<=(theta_bins-1) && ir_max>=0 && ir_max<=(rho_bins-1)) {
		 
		 for(int ix=0; ix<3; ix++){
			for(int iy=0; iy<3; iy++){
			
				if((it_max+ix-1)<0 || (it_max+ix-1)>(theta_bins-1)) continue;
				if((ir_max+iy-1)<0 || (ir_max+iy-1)>(rho_bins-1)) continue;
				
				if(!(ix==1&&iy==1)){
				for(int ic=0; ic<boxes[it_max+ix-1][ir_max+iy-1]; ic++){
				
					fitparam params1;
				
					params1.x1 = x1_val[it_max+ix-1][ir_max+iy-1][ic];
					params1.x2 = x2_val[it_max+ix-1][ir_max+iy-1][ic];
					params1.z1 = z1_val[it_max+ix-1][ir_max+iy-1][ic];
					params1.z2 = z2_val[it_max+ix-1][ir_max+iy-1][ic];
					params1.theta = theta_val[it_max+ix-1][ir_max+iy-1][ic];
					params1.rho = rho_val[it_max+ix-1][ir_max+iy-1][ic];
					
					final_xpars.push_back(params1);
					
					}
				}
			
				deno += boxes[it_max+ix-1][ir_max+iy-1];
				t_s += (theta_min+theta_binsize*(it_max+ix-1)) * boxes[it_max+ix-1][ir_max+iy-1];
				r_s += (rho_min+rho_binsize*(ir_max+iy-1)) * boxes[it_max+ix-1][ir_max+iy-1];
				
				}
			}		
		}
		
		 t_s *= 1./deno;
		 r_s *= 1./deno;
		 
	}
	else{ break; }
 }
 
 // optimize in y-z plane
 
	for(int iloop=0; iloop<nloop; iloop++){
	 
		if(fabs(ty_s - ty_opt) > 2*theta_binsize || fabs(ry_s - ry_opt) > 2*rho_binsize){
		 
		 if(iloop > 0){
		 
			for(int it_1=0; it_1<theta_bins; it_1++){
				for(int jc_1=0; jc_1<rho_bins; jc_1++){
		 
						if( (ty_s > (theta_min + (it_1)*theta_binsize)) && (ty_s < (theta_min + (it_1+1)*theta_binsize))
							&& (ry_s > (rho_min + (jc_1)*rho_binsize)) && (ry_s < (rho_min + (jc_1+1)*rho_binsize)) )
						{	
							ity_max = it_1; iry_max = jc_1;
						}
					}
				}
			
			ty_opt = ty_s;
			ry_opt = ry_s;
		
			}
			
		 deno1 = 0;
		 ty_s =0; ry_s = 0;
		 
		 if(ity_max>=0 && ity_max<=(theta_bins-1) && iry_max>=0 && iry_max<=(rho_bins-1)) {
		 
		 for(int ix=0; ix<3; ix++){
			for(int iy=0; iy<3; iy++){
			
				if((ity_max+ix-1)<0 || (ity_max+ix-1)>(theta_bins-1)) continue;
				if((iry_max+iy-1)<0 || (iry_max+iy-1)>(rho_bins-1)) continue;
			
				deno1 += boyes[ity_max+ix-1][iry_max+iy-1];
				ty_s += (theta_min+theta_binsize*(ity_max+ix-1)) * boyes[ity_max+ix-1][iry_max+iy-1];
				ry_s += (rho_min+rho_binsize*(iry_max+iy-1)) * boyes[ity_max+ix-1][iry_max+iy-1];
				
				if(!(ix==1&&iy==1)){
				for(int ic=0; ic<boyes[ity_max+ix-1][iry_max+iy-1]; ic++){
			
					fitparam params1;
					
					params1.x1 = y1_val[ity_max+ix-1][iry_max+iy-1][ic];
					params1.x2 = y2_val[ity_max+ix-1][iry_max+iy-1][ic];
					params1.z1 = zy1_val[ity_max+ix-1][iry_max+iy-1][ic];
					params1.z2 = zy2_val[ity_max+ix-1][iry_max+iy-1][ic];
					params1.theta = theta_yval[ity_max+ix-1][iry_max+iy-1][ic];
					params1.rho = rho_yval[ity_max+ix-1][iry_max+iy-1][ic];
					
					final_ypars.push_back(params1);
					
					}
				}
			
			}
		}		
	}
	
	ty_s *= 1./deno1;
	ry_s *= 1./deno1;
		 
	}
	else{ break; }
 }
	
#endif


	r_s = (fabs(sin(t_s)) > 1.e-2)?(r_s*1./sin(t_s)):1.e-2;
	ry_s = (fabs(sin(ty_s)) > 1.e-2)?(ry_s*1./sin(ty_s)):1.e-2;
	
	
	//  Find line points & line in x-z plane //

	vector<point> final_xpoints;
	
	// cout<<"Number of pairs "<<final_xpars.size()<<endl;
	for(int ic=0; ic<final_xpars.size(); ic++){

		point point_a, point_b;
		point_a.x = final_xpars[ic].x1;
		point_a.z = final_xpars[ic].z1;
		point_b.x = final_xpars[ic].x2;
		point_b.z = final_xpars[ic].z2;
		
		bool p1_pass = true;
		if(final_xpoints.size()>0){
		for(int jd=0; jd<final_xpoints.size(); jd++){
				if(fabs(final_xpoints[jd].x - point_a.x) < 1.e-6 &&  fabs(final_xpoints[jd].z - point_a.z) < 1.e-6){
					p1_pass = false;
					break;
					}
			}
		}
		if(p1_pass) { final_xpoints.push_back(point_a); }
		
		bool p2_pass = true;
		if(final_xpoints.size()>0){
		for(int jd=0; jd<final_xpoints.size(); jd++){
				if(fabs(final_xpoints[jd].x - point_b.x) < 1.e-6 &&  fabs(final_xpoints[jd].z - point_b.z) < 1.e-6){
					p2_pass = false;
					break;
					}
			}
		}
		if(p2_pass) { final_xpoints.push_back(point_b); }
		
	}
	
	if(final_xpoints.size() < 2) return;

	double xzsum =0;
	double xsum = 0;
	double x2sum = 0;
	double zsum = 0;
	
	for(int ip=0; ip<final_xpoints.size(); ip++){
//		// cout<<"point:"<<ip+1<<" x:"<<final_xpoints[ip].x<<" z:"<<final_xpoints[ip].z<<endl;
		xzsum += (final_xpoints[ip].x * final_xpoints[ip].z);
		xsum += final_xpoints[ip].x;
		x2sum += (final_xpoints[ip].x * final_xpoints[ip].x);
		zsum += final_xpoints[ip].z;
		}
	
//	double c_final = 	(zsum*xsum - x2sum*zsum)*1./(xsum*xsum - final_xpoints.size()*x2sum);
//	double m_final = (zsum - c_final*final_xpoints.size())*1./xsum;
	
	double m_initial = (final_xpoints.size()*xzsum - xsum*zsum)*1./max(final_xpoints.size()*x2sum - xsum*xsum,1.e-3);
	double c_initial = (zsum - m_initial*xsum)*1./final_xpoints.size();
	
	// cout<<"m_inital "<<m_initial<<" c_initial "<<c_initial<<endl;
	
	double chi2 = 0;
	for(int ip=0; ip<final_xpoints.size(); ip++){
		
		double ypred = (m_initial*final_xpoints[ip].x+c_initial);
		double yerr = m_initial*strip_width*1./sqrt(12);
		double chi2_try = (final_xpoints[ip].z - ypred)*(final_xpoints[ip].z - ypred)*1./(yerr*yerr);
		
		chi2 += chi2_try;
		}	
	
	// cout<<"chi2/NDF:"<<chi2<<"/"<<final_xpoints.size()<<endl;

	for(int ip=0; ip<nplanes; ip++){
	
		double tryx = (zval[ip] - c_initial)*1./m_initial;
		double min_dis = 1000;
		int npoints_plane = 0;
	
		for(int fp=0; fp<final_xpoints.size(); fp++){
			
			if(fabs(zval[ip] - final_xpoints[fp].z) < 1.e-6){
				
				double xdiff =  fabs(tryx - final_xpoints[fp].x);
				
				if(xdiff < min_dis) { min_dis = xdiff; }
				npoints_plane++;
			}
		}
		
		if(npoints_plane>1){
			for(int fp=0; fp<final_xpoints.size(); fp++){
				if(fabs(zval[ip] - final_xpoints[fp].z) < 1.e-6){
					double xdiff =  fabs(tryx - final_xpoints[fp].x);
					if(xdiff > min_dis){
						final_xpoints.erase(final_xpoints.begin()+fp);
						}
					}
				}
			}
		
	}
	
  double m_final, c_final;	
	
  for(int iloop=0; iloop<10; iloop++){
	
	for(int ip=0; ip<nplanes; ip++){
		
		double tryx = (zval[ip] - c_initial)*1./m_initial;
		
		for(int fp=0; fp<final_xpoints.size(); fp++){
			
			if(fabs(zval[ip] - final_xpoints[fp].z) < 1.e-6){
				
				double xdiff =  fabs(tryx - final_xpoints[fp].x);
				if(fabs(xdiff) > 10){
					final_xpoints.erase(final_xpoints.begin()+fp);
				}
			}
		}	
	}
	
	xzsum = 0; xsum = 0; x2sum = 0; zsum = 0;
	
	for(int ip=0; ip<final_xpoints.size(); ip++){
//		// cout<<"point:"<<ip+1<<" x:"<<final_xpoints[ip].x<<" z:"<<final_xpoints[ip].z<<endl;
		xzsum += (final_xpoints[ip].x * final_xpoints[ip].z);
		xsum += final_xpoints[ip].x;
		x2sum += (final_xpoints[ip].x * final_xpoints[ip].x);
		zsum += final_xpoints[ip].z;
		}
	
	m_final = (final_xpoints.size()*xzsum - xsum*zsum)*1./max(final_xpoints.size()*x2sum - xsum*xsum,1.e-3);
	c_final = (zsum - m_final*xsum)*1./final_xpoints.size();
	
	// cout<<"m_final "<<m_final<<" c_final "<<c_final<<endl;
	
	if((fabs(m_final - m_initial) < 0.2*fabs(m_final+0.001)) && (fabs(c_final - c_initial) < 0.2*fabs(c_final+0.001))) break;
	
	m_initial = m_final;
	c_initial = c_final;
	
	}
	
	chi2 = 0;
	for(int ip=0; ip<final_xpoints.size(); ip++){
		
		double ypred = (m_final*final_xpoints[ip].x+c_final);
		double yerr = m_final*strip_width*1./sqrt(12);
		double chi2_try = (final_xpoints[ip].z - ypred)*(final_xpoints[ip].z - ypred)*1./(yerr*yerr);
		
		chi2 += chi2_try;
		}	
	
	// cout<<"chi2/NDF:"<<chi2<<"/"<<final_xpoints.size()-2<<endl;
	
	 cout<<"Number of points "<<final_xpoints.size()<<endl;
//	 cout<<"And points are: "<<endl;
	for(int ip=0; ip<final_xpoints.size(); ip++){
//		cout<<"point:"<<ip+1<<" x:"<<final_xpoints[ip].x<<" z:"<<final_xpoints[ip].z<<endl;
	}
	
	
	// x-z plane ends //
	
	// find points and line in y-z plane
	
	// cout<<"\n Results: y-z \n"<<endl;
	
	vector<point> final_ypoints;
	
	// cout<<"Number of pairs "<<final_ypars.size()<<endl;
	if(final_ypars.size()<1) return;
	
	for(int ic=0; ic<final_ypars.size(); ic++){

		point point_a, point_b;
		point_a.x = final_ypars[ic].x1;
		point_a.z = final_ypars[ic].z1;
		point_b.x = final_ypars[ic].x2;
		point_b.z = final_ypars[ic].z2;
		
		bool p1_pass = true;
		if(final_ypoints.size()>0){
		for(int jd=0; jd<final_ypoints.size(); jd++){
				if(fabs(final_ypoints[jd].x - point_a.x) < 1.e-6 &&  fabs(final_ypoints[jd].z - point_a.z) < 1.e-6){
					p1_pass = false;
					break;
					}
			}
		}
		if(p1_pass) { final_ypoints.push_back(point_a); }
		
		bool p2_pass = true;
		if(final_ypoints.size()>0){
		for(int jd=0; jd<final_ypoints.size(); jd++){
				if(fabs(final_ypoints[jd].x - point_b.x) < 1.e-6 &&  fabs(final_ypoints[jd].z - point_b.z) < 1.e-6){
					p2_pass = false;
					break;
					}
			}
		}
		if(p2_pass) { final_ypoints.push_back(point_b); }
		
	}
	
	if(final_ypoints.size() < 2) return;

	double yzsum =0;
	double z1sum = 0;
	double y2sum = 0;
	double ysum = 0;
	
	for(int ip=0; ip<final_ypoints.size(); ip++){
		yzsum += (final_ypoints[ip].x * final_ypoints[ip].z);
		ysum += final_ypoints[ip].x;
		y2sum += (final_ypoints[ip].x * final_ypoints[ip].x);
		z1sum += final_ypoints[ip].z;
		}
	
	double my_initial = (final_ypoints.size()*yzsum - ysum*z1sum)*1./max(final_ypoints.size()*y2sum - ysum*ysum,1.e-3);
	double cy_initial = (z1sum - my_initial*ysum)*1./final_ypoints.size();
	
	// cout<<"my_inital "<<my_initial<<" cy_initial "<<cy_initial<<endl;
	
	double ychi2 = 0;
	for(int ip=0; ip<final_ypoints.size(); ip++){
		
		double ypred = (my_initial*final_ypoints[ip].x+cy_initial);
		double yerr = my_initial*strip_width*1./sqrt(12);
		double chi2_try = (final_ypoints[ip].z - ypred)*(final_ypoints[ip].z - ypred)*1./(yerr*yerr);
		
		ychi2 += chi2_try;
		}	
	
	
//	// cout<<"chi2/NDF:"<<chi2<<"/"<<final_ypoints.size()<<endl;

	for(int ip=0; ip<nplanes; ip++){
	
		double tryy = (zval[ip] - cy_initial)*1./my_initial;
		double min_dis = 1000;
		int npoints_plane = 0;
	
		for(int fp=0; fp<final_ypoints.size(); fp++){
			
			if(fabs(zval[ip] - final_ypoints[fp].z) < 1.e-6){
				
				double ydiff =  fabs(tryy - final_ypoints[fp].x);
				
				if(ydiff < min_dis) { min_dis = ydiff; }
				npoints_plane++;
			}
		}
		
		if(npoints_plane>1){
			for(int fp=0; fp<final_ypoints.size(); fp++){
				if(fabs(zval[ip] - final_ypoints[fp].z) < 1.e-6){
					double ydiff =  fabs(tryy - final_ypoints[fp].x);
					if(ydiff > min_dis){
						final_ypoints.erase(final_ypoints.begin()+fp);
						}
					}
				}
			}
		
	}
	
  double my_final, cy_final;	
	
  for(int iloop=0; iloop<10; iloop++){
	
	for(int ip=0; ip<nplanes; ip++){
		
		double tryy = (zval[ip] - cy_initial)*1./my_initial;
		
		for(int fp=0; fp<final_ypoints.size(); fp++){
			
			if(fabs(zval[ip] - final_ypoints[fp].z) < 1.e-6){
				
				double ydiff =  fabs(tryy - final_ypoints[fp].x);
				if(fabs(ydiff) > 10){
					final_ypoints.erase(final_ypoints.begin()+fp);
				}
			}
		}	
	}
	
	yzsum = 0; ysum = 0; y2sum = 0; z1sum = 0;
	
	for(int ip=0; ip<final_ypoints.size(); ip++){
		yzsum += (final_ypoints[ip].x * final_ypoints[ip].z);
		ysum += final_ypoints[ip].x;
		y2sum += (final_ypoints[ip].x * final_ypoints[ip].x);
		z1sum += final_ypoints[ip].z;
		}
	
	my_final = (final_ypoints.size()*yzsum - ysum*z1sum)*1./max(final_ypoints.size()*y2sum - ysum*ysum,1.e-3);
	cy_final = (z1sum - my_final*ysum)*1./final_ypoints.size();
	
	// cout<<"my_final "<<my_final<<" cy_final "<<cy_final<<endl;
	
	if((fabs(my_final - my_initial) < 0.2*fabs(my_final+0.001)) && (fabs(cy_final - cy_initial) < 0.2*fabs(cy_final+0.001))) break;
	
	my_initial = my_final;
	cy_initial = cy_final;
	
	}
	
	ychi2 = 0;
	for(int ip=0; ip<final_ypoints.size(); ip++){
			
		double ypred = (my_final*final_ypoints[ip].x+cy_final);
		double yerr = my_final*strip_width*1./sqrt(12);
		double chi2_try = (final_ypoints[ip].z - ypred)*(final_ypoints[ip].z - ypred)*1./(yerr*yerr);
	
		ychi2 += chi2_try; 
		}	
	
	// cout<<"chi2/NDF:"<<ychi2<<"/"<<final_ypoints.size()-2<<endl;
	
	// cout<<"Number of points "<<final_ypoints.size()<<endl;
	// cout<<"And points are: "<<endl;
	for(int ip=0; ip<final_ypoints.size(); ip++){
		// cout<<"point:"<<ip+1<<" y:"<<final_ypoints[ip].x<<" z:"<<final_ypoints[ip].z<<endl;
	}
	
	
	// y-z plane ends
	*/
}
