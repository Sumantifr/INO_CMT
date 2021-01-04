#include "Cosmic_trig.h"

/*
#define Remove_Noise
//#define PRINT_Details
#define Optimize
#define Optimize_Loop
#define x_z_only
*/

using namespace std;

//bool Cosmic_trig(bool **xhitx_strips, bool **yhitx_strips)
void Cosmic_trig_top(const bool xhitx_strips[nstrips][nplanes],const bool yhitx_strips[nstrips][nplanes], const float zval[nplanes], bool &xzpass)
{
	
	bool xz_pass = false;
	
	float xval[nstrips][nplanes];
	float yval[nstrips][nplanes];
//	float zval[nplanes];
	
	for(int ip=0; ip<nplanes; ip++){
//		zval[ip] = plane_gap*(ip);
		
		for(int js=0; js<nstrips; js++){
//			#pragma HLS pipeline II=1
			xval[js][ip] = (js)*strip_gap + 0.5*(1+2*js)*strip_width - (0.5*nstrips*strip_width) - (0.5*nstrips-1)*strip_gap - 0.5*strip_gap;
			yval[js][ip] = (js)*strip_gap + 0.5*(1+2*js)*strip_width - (0.5*nstrips*strip_width) - (0.5*nstrips-1)*strip_gap - 0.5*strip_gap;
			}
		
		}

	#ifdef Remove_Noise
	
	// remove noise in x strips
	
	for(int ip=0; ip<nplanes; ip++){
		#pragma HLS unroll factor=8
		// if there are 4 strips with hits in a single layer, discard the layer //
		
		int nstr_hit = 0; 
		for(int ix1=0; ix1<nstrips; ix1++){
			if(xhitx_strips[ix1][ip]){
				nstr_hit++;
			}
		}
		if(nstr_hit>=5){
			for(int ix1=0; ix1<nstrips; ix1++){
				xhitx_strips[ix1][ip] = 0;
				}
			}
		
		
		// done//
		
		// if there are hits in 3 consecutive strips, merge them //
		
		for(int ix1=0; ix1<(nstrips-2); ix1++){
			if(xhitx_strips[ix1][ip] && xhitx_strips[ix1+1][ip] && xhitx_strips[ix1+2][ip]){
				xhitx_strips[ix1][ip] = xhitx_strips[ix1+2][ip] = 0;
				}
		}
		
		//  if there are hits in 2 consecutive strips, merge them //
		
		for(int ix1=0; ix1<(nstrips-1); ix1++){
			if(xhitx_strips[ix1][ip] && xhitx_strips[ix1+1][ip]){
				xhitx_strips[ix1+1][ip] = 0;
				xval[ix1][ip] = 0.5*(xval[ix1][ip]+xval[ix1+1][ip]);
				}
		}
		
		// done
	}//
	
	#endif
	
	
	const int theta_bins = 10;//40;
//	float theta_max = M_PI + 0.05; float theta_min = 0 - 0.05;
	float theta_max = 0.5*M_PI + 0.05; float theta_min = -0.5*M_PI;
	float theta_binsize = (theta_max-theta_min)*1./theta_bins;
	
	const int rho_bins = 15;//30;
	float rho_max = 28.+ 1; float rho_min = -28-1;
	float rho_binsize = (rho_max-rho_min)*1./rho_bins;
	
	static const int max_entry = 100;
	
	int boxes[theta_bins][rho_bins] = {{0}};
	
	int npf_points = 0;
	point prefinal[theta_bins][rho_bins][max_entry];

	//// TRACK Finder (Method: Hough transformation)  /////
	
	// fill grids in x-z plane   
	
	for(int ip=0; ip<nplanes; ip++){
		
		for(int ix1=0; ix1<nstrips; ix1++){
			#pragma HLS unroll factor=4
		//	#pragma HLS pipeline II=1
			
			if(xhitx_strips[ix1][ip]){
				
				for(int jp=(ip+1); jp<nplanes; jp++){
					
					for(int ix2=0; ix2<nstrips; ix2++){
						#pragma HLS pipeline II=1
						
						if(xhitx_strips[ix2][jp]){
							
							float theta_test, tantheta_test, rho_test;
							
							tantheta_test = -1.*(xval[ix2][jp] - xval[ix1][ip])*1./(zval[jp] - zval[ip]);
							theta_test = atan(tantheta_test);
							rho_test = (xval[ix1][ip] + tantheta_test * zval[ip])*1./sqrt(1+tantheta_test*tantheta_test);
							rho_test += (xval[ix2][jp] + tantheta_test * zval[jp])*1./sqrt(1+tantheta_test*tantheta_test);
							rho_test = 0.5* rho_test;
						
							if(rho_test>=rho_min && rho_test<rho_max && theta_test>=theta_min && theta_test<theta_max){
								
								int itheta_bin = int((theta_test-theta_min)*1/theta_binsize);
								int irho_bin =	int((rho_test-rho_min)*1./rho_binsize);
								boxes[itheta_bin][irho_bin] += 1;
								
								point pa;
								pa.x = xval[ix1][ip]; pa.z = zval[ip];
								prefinal[itheta_bin][irho_bin][2*(boxes[itheta_bin][irho_bin]-1)] = pa;
								pa.x = xval[ix2][jp]; pa.z = zval[jp];
								prefinal[itheta_bin][irho_bin][2*(boxes[itheta_bin][irho_bin]-1)+1] = pa;
								
							}
							
							
						}
					
					}
				}
				
			}
			
		}
	}
	
	
	// find the grid with the highest population

	int itx_max = -1; int irx_max = -1; 
	
	float max_xcon = -100;

	for(int it=0; it<theta_bins; it++){
		for(int jc=0; jc<rho_bins; jc++){
			#pragma HLS pipeline II=1
			if(boxes[it][jc]>= max_xcon){
				max_xcon = boxes[it][jc];
				itx_max = it;
				irx_max = jc;
				}
			}
		}
	
	// find pairs belonging to the grid with the highest population

	float deno, tx_s, rx_s;
	
	tx_s = -100; rx_s =  -100;
	rx_s = (fabs(sin(tx_s)) > 1.e-2)?(rx_s*1./sin(tx_s)):1.e-2;
	
//	cout<<"==== x-z ==== "<<endl;

   //  Find points to be fitted by a line in x-z plane //

	point final_xpoints[20];
	int nxpoint = 0;
	
	int npfxpoint = 2*boxes[itx_max][irx_max];
	
	// option 1
	
	for(int ic=0; ic<(npfxpoint); ic++){
		#pragma HLS loop_tripcount min=1 max=100
		 
		 if(nxpoint==0){
			final_xpoints[nxpoint] = prefinal[itx_max][irx_max][ic];
			nxpoint++;
			}
		
		 if(nxpoint>0){
			 bool match_found = false;
			 for(int jc=0; jc<(nxpoint); jc++){
				#pragma HLS loop_tripcount min=1 max=100
				if( (fabs(prefinal[itx_max][irx_max][ic].x - final_xpoints[jc].x) < 1.e-6) && (fabs(prefinal[itx_max][irx_max][ic].z - final_xpoints[jc].z) < 1.e-6) ){
					match_found = true;
					break;
					}	  
				 }
				if(!match_found){
					 final_xpoints[nxpoint] = prefinal[itx_max][irx_max][ic];
					 nxpoint++;
					} 
			  }
		 
	}
	
	/*
	nxpoint = 7;
	final_xpoints[0].x = -10.5;	final_xpoints[0].z = 0;	
	final_xpoints[1].x = -10.5;	final_xpoints[1].z = 4;	
	final_xpoints[2].x = -10.5;	final_xpoints[2].z = 8;	
	final_xpoints[3].x = -10.5;	final_xpoints[3].z = 12;	
	final_xpoints[4].x = -7.5;	final_xpoints[4].z = 12;	
	final_xpoints[5].x = -7.5;	final_xpoints[5].z = 16;	
	final_xpoints[6].x = -7.5;	final_xpoints[6].z = 20;	
	*/

	// option 1 ends
	
	// option 2
	 /*
	for(int ic=0; ic<(npfxpoint); ic++){
		#pragma HLS loop_tripcount min=1 max=100
		final_xpoints[nxpoint] = prefinal[itx_max][irx_max][ic];
		nxpoint++;
	}
	 
	for(int ic=0; ic<(nxpoint); ic++){
		#pragma HLS loop_tripcount min=1 max=100
		for(int jc=(ic+1); jc<(nxpoint); jc++){
			#pragma HLS pipeline II=1
			#pragma HLS loop_tripcount min=1 max=100
			if( (fabs(final_xpoints[ic].x - final_xpoints[jc].x) < 1.e-6) && (fabs(final_xpoints[ic].z - final_xpoints[jc].z) < 1.e-6) ){
				for(int kd=jc; kd<(nxpoint-1); kd++){
					#pragma HLS loop_tripcount min=1 max=50
					final_xpoints[kd] = final_xpoints[kd+1];
					}
				jc -= 1;
				nxpoint--;
				}
		}
	}
	*/
	// option 2 ends
	
	
	
	if(nxpoint < 2) return;// xzpass;
	
	//// TRACK Finder ends  /////

	//// Track fitter starts ///

	float xzsum =0;
	float xsum = 0;
	float x2sum = 0;
	float zsum = 0;
	
	// initial fit //
	
	for(int ip=0; ip<nxpoint; ip++){
		#pragma HLS loop_tripcount min=1 max=28
//		#pragma HLS pipeline II=1
		xzsum += (final_xpoints[ip].x * final_xpoints[ip].z);
		xsum += final_xpoints[ip].x;
		x2sum += (final_xpoints[ip].x * final_xpoints[ip].x);
		zsum += final_xpoints[ip].z;
		}
	
	float m_initial = (nxpoint*xzsum - xsum*zsum)*1./max(nxpoint*x2sum - xsum*xsum,float(1.e-3));
	float c_initial = (zsum - m_initial*xsum)*1./nxpoint;
	
	// initial chi2 //
	
	float chi2 = 0;
	for(int ip=0; ip<nxpoint; ip++){
		#pragma HLS loop_tripcount min=1 max=28
//		#pragma HLS pipeline II=1
		float ypred = (m_initial*final_xpoints[ip].x + c_initial);
		float yerr = m_initial*strip_width*1./sqrt(12);
		float chi2_try = (final_xpoints[ip].z - ypred)*(final_xpoints[ip].z - ypred)*1./(yerr*yerr);
		
		chi2 += chi2_try;
		}	
	
   // remove outliers //

	for(int ip=0; ip<nplanes; ip++){
		
		float tryx = (zval[ip] - c_initial)*1./m_initial;
		float min_dis = 10000;
		int npoints_plane = 0;
	
		for(int fp=0; fp<nxpoint; fp++){
//			#pragma HLS pipeline II=1
			#pragma HLS loop_tripcount min=1 max=28
			
			if(fabs(zval[ip] - final_xpoints[fp].z) < 1.e-6){
				
				float xdiff =  fabs(tryx - final_xpoints[fp].x);
				if(xdiff < min_dis) { min_dis = xdiff; }
				npoints_plane++;
			}
		
		}
		
		if(npoints_plane>1){
			
			int which_points[nstrips];
			int nextra = 0;
			
			for(int fp=0; fp<nxpoint; fp++){
				#pragma HLS loop_tripcount min=1 max=28
				if(fabs(zval[ip] - final_xpoints[fp].z) < 1.e-6){
					float xdiff =  fabs(tryx - final_xpoints[fp].x);
					if(xdiff > min_dis){
						which_points[nextra] = fp;
						nextra++;
						}
					}
				}
			
			for(int ip=0; ip<nextra; ip++){
				#pragma HLS loop_tripcount min=1 max=28
				for(int fp=which_points[ip]; fp<nxpoint; fp++){
					#pragma HLS loop_tripcount min=1 max=28
						final_xpoints[fp] = final_xpoints[fp+1];
					}
				nxpoint -= 1;
				which_points[ip+1] -= 1;
				}
			
			}
		
	}
	 // remove outliers done //

 // final fit //

	
  float m_final, c_final;	
  int niter = 0;	
  if(niter>0){
     // switched off to save time
	for(int iloop=0; iloop<niter; iloop++){
	
	// removing points too far from line
	
		for(int ip=0; ip<nplanes; ip++){
//			#pragma HLS pipeline II=1
			float tryx = (zval[ip] - c_initial)*1./m_initial;
		
			int nextra=0;
			int which_points[nstrips];
		
			for(int fp=0; fp<nxpoint; fp++){
			
				#pragma HLS loop_tripcount min=1 max=28
				if(fabs(zval[ip] - final_xpoints[fp].z) < 1.e-6){
				
					float xdiff =  fabs(tryx - final_xpoints[fp].x);
					if(fabs(xdiff) > 10){
						which_points[nextra] = fp;
						nextra++;
					}
				}
			}
		
			for(int ip=0; ip<nextra; ip++){
				#pragma HLS loop_tripcount min=1 max=28
				for(int fp=which_points[ip]; fp<nxpoint; fp++){
					#pragma HLS loop_tripcount min=1 max=28
					final_xpoints[fp] = final_xpoints[fp+1];
					}
				nxpoint -= 1;
				which_points[ip+1] -= 1;
			}	
		}
	
		xzsum = 0; xsum = 0; x2sum = 0; zsum = 0;
	
		for(int ip=0; ip<nxpoint; ip++){
			#pragma HLS loop_tripcount min=1 max=28
//			#pragma HLS unroll factor=2
			xzsum += (final_xpoints[ip].x * final_xpoints[ip].z);
			xsum += final_xpoints[ip].x;
			x2sum += (final_xpoints[ip].x * final_xpoints[ip].x);
			zsum += final_xpoints[ip].z;
			}
	
		m_final = (nxpoint*xzsum - xsum*zsum)*1./max(nxpoint*x2sum - xsum*xsum,float(1.e-3));
		c_final = (zsum - m_final*xsum)*1./nxpoint;
	
		if((fabs(m_final - m_initial) < 0.2*fabs(m_final+0.001)) && (fabs(c_final - c_initial) < 0.2*fabs(c_final+0.001))) break;
	
		m_initial = m_final;
		c_initial = c_final;
	
		}
  
   }
	
	if(niter==0){
		
		xzsum = 0; xsum = 0; x2sum = 0; zsum = 0;
	
		for(int ip=0; ip<nxpoint; ip++){
		#pragma HLS loop_tripcount min=1 max=28
//		#pragma HLS pipeline II=1
			xzsum += (final_xpoints[ip].x * final_xpoints[ip].z);
			xsum += final_xpoints[ip].x;
			x2sum += (final_xpoints[ip].x * final_xpoints[ip].x);
			zsum += final_xpoints[ip].z;
		}
	
		m_final = (nxpoint*xzsum - xsum*zsum)*1./max(nxpoint*x2sum - xsum*xsum,float(1.e-3));
		c_final = (zsum - m_final*xsum)*1./nxpoint;
		
		}
	
	// final chi^2
	
	chi2 = 0;
	for(int ip=0; ip<nxpoint; ip++){
		#pragma HLS loop_tripcount min=1 max=28
//		#pragma HLS pipeline II=1
		float ypred = (m_final*final_xpoints[ip].x + c_final);
		float yerr = m_final*strip_width*1./sqrt(12);
		float chi2_try = (final_xpoints[ip].z - ypred)*(final_xpoints[ip].z - ypred)*1./(yerr*yerr);
		
		chi2 += chi2_try;
		}	
	
	 cout<<"Final chi2/NDF:"<<chi2<<"/"<<nxpoint-2<<endl;
	
	 cout<<"Number of points "<<nxpoint<<endl;
	 cout<<"And points are: "<<endl;
	 for(int ip=0; ip<nxpoint; ip++){
		#pragma HLS loop_tripcount min=1 max=28
		cout<<"point:"<<ip+1<<" x:"<<final_xpoints[ip].x<<" z:"<<final_xpoints[ip].z<<endl;
	 }
	
	float chi2x = chi2*1./(nxpoint-2);
	
	if(nxpoint >= 4 && chi2x<5.) {  
		xz_pass = true; 
		}

	xzpass = xz_pass;

	// x-z plane ends //
	
	// y-z plane starts //
	
	/*
	 
	#ifdef Remove_Noise 
	  
	// remove noise in y-strips
	
	for(int ip=0; ip<nplanes; ip++){
		#pragma HLS unroll factor=1
		// if there are 4 strips with hits in a single layer, discard the layer //
		
		int nstr_hit = 0; 
		
		for(int iy1=0; iy1<nstrips; iy1++){
			if(yhitx_strips[iy1][ip]){
				nstr_hit++;
			}
		}
		if(nstr_hit>=5){
			for(int iy1=0; iy1<nstrips; iy1++){
				yhitx_strips[iy1][ip] = 0;
				}
			}
		
			
		// done//
		
		// if there are hits in 3 consecutive strips, merge them //
		
		for(int iy1=0; iy1<(nstrips-2); iy1++){
			if(yhitx_strips[iy1][ip] && yhitx_strips[iy1+1][ip] && yhitx_strips[iy1+2][ip]){
				yhitx_strips[iy1][ip] = yhitx_strips[iy1+2][ip] = 0;
				}
		}
		
		//  if there are hits in 2 consecutive strips, merge them //
		
		for(int iy1=0; iy1<(nstrips-1); iy1++){
			if(yhitx_strips[iy1][ip] && yhitx_strips[iy1+1][ip]){
				yhitx_strips[iy1+1][ip] = 0;
				yval[iy1][ip] = 0.5*(yval[iy1][ip]+yval[iy1+1][ip]);
				}
		}
		
		// done
	}//
	
	#endif
	
	float theta_yval[theta_bins][rho_bins][max_entry];
	float rho_yval[theta_bins][rho_bins][max_entry];
	
	float y1_val[theta_bins][rho_bins][max_entry];
	float y2_val[theta_bins][rho_bins][max_entry];
	float zy1_val[theta_bins][rho_bins][max_entry];
	float zy2_val[theta_bins][rho_bins][max_entry];
	
	
	int boyes[theta_bins][rho_bins] = {{0}};
	
	
	cout<<"==== y-z ==== "<<endl;
	
	// fill grids in y-z plane
	
	for(int ip=0; ip<nplanes; ip++){
		
		for(int iy1=0; iy1<nstrips; iy1++){
			
			if(yhitx_strips[iy1][ip]){
				
				for(int jp=(ip+1); jp<nplanes; jp++){
				#pragma HLS loop_tripcount min=1
//					if(jp>(ip+3)) continue; // forming pairs only between 3 consecutive layers
				
					for(int iy2=0; iy2<nstrips; iy2++){
					
						if(yhitx_strips[iy2][jp]){
							
							float theta_test, tantheta_test, rho_test;
							
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
								theta_yval[itheta_bin][irho_bin][boxes[itheta_bin][irho_bin]-1] = theta_test;
								rho_yval[itheta_bin][irho_bin][boxes[itheta_bin][irho_bin]-1] = rho_test;
								y1_val[itheta_bin][irho_bin][boxes[itheta_bin][irho_bin]-1] = yval[iy1][ip];
								y2_val[itheta_bin][irho_bin][boxes[itheta_bin][irho_bin]-1] = yval[iy2][jp];
								zy1_val[itheta_bin][irho_bin][boxes[itheta_bin][irho_bin]-1] = zval[ip];
								zy2_val[itheta_bin][irho_bin][boxes[itheta_bin][irho_bin]-1] = zval[jp];
								
							}
							
						}
					
					}
				}
				
			}
			
		}
	}
	
	// find the grid with the highest population
	
	int ity_max = -1; int iry_max = -1;
	
	float max_ycon = -100;
	
	for(int it=0; it<theta_bins; it++){
		for(int jc=0; jc<rho_bins; jc++){
			if(boyes[it][jc]>= max_ycon){
				max_ycon = boyes[it][jc];
				ity_max = it;
				iry_max = jc;
				}
			}
		}	
	
	// pairs belonging to the grid with the highest population
	
	fitparam paramsy[20];
	int npairsy = boyes[ity_max][iry_max];
		
	for(int ic=0; ic<boyes[ity_max][iry_max]; ic++){
		#pragma HLS loop_tripcount min=1 max=100
		paramsy[ic].x1 = y1_val[ity_max][iry_max][ic];
		paramsy[ic].x2 = y2_val[ity_max][iry_max][ic];
		paramsy[ic].z1 = zy1_val[ity_max][iry_max][ic];
		paramsy[ic].z2 = zy2_val[ity_max][iry_max][ic];
		paramsy[ic].theta = theta_yval[ity_max][iry_max][ic];
		paramsy[ic].rho = rho_yval[ity_max][iry_max][ic];
		
		}	
	
	float ty_opt = 	(theta_min+theta_binsize*(ity_max+0.5));
	float ry_opt = 	(rho_min+rho_binsize*(iry_max+0.5));
	
	float deno1, ty_s, ry_s;
	
	ty_s = -100; ry_s =  -100;
	
	ry_s = (fabs(sin(ty_s)) > 1.e-2)?(ry_s*1./sin(ty_s)):1.e-2;
	
	// find points and line in y-z plane
	
	point final_ypoints[20];
	int nypoint = 0;
	
	// cout<<"Number of pairs "<<final_ypars.size()<<endl;
	if(npairsy<1) return;
	
	for(int ic=0; ic<npairsy; ic++){
		#pragma HLS loop_tripcount min=1 max=28
		
		point point_a, point_b;
		point_a.x = paramsy[ic].x1;
		point_a.z = paramsy[ic].z1;
		point_b.x = paramsy[ic].x2;
		point_b.z = paramsy[ic].z2;
		
		bool p1_pass = true;
		if(nypoint>0){
		for(int jd=0; jd<nypoint; jd++){
			#pragma HLS loop_tripcount min=1 max=28
				if(fabs(final_ypoints[jd].x - point_a.x) < 1.e-6 &&  fabs(final_ypoints[jd].z - point_a.z) < 1.e-6){
					p1_pass = false;
					break;
					}
			}
		}
		if(p1_pass) { final_ypoints[nypoint] = point_a; nypoint++;}
		
		bool p2_pass = true;
		if(nypoint>0){
		for(int jd=0; jd<nypoint; jd++){
			#pragma HLS loop_tripcount min=1 max=28
				if(fabs(final_ypoints[jd].x - point_b.x) < 1.e-6 &&  fabs(final_ypoints[jd].z - point_b.z) < 1.e-6){
					p2_pass = false;
					break;
					}
			}
		}
		if(p2_pass) { final_ypoints[nypoint] = point_b; nypoint++;}
		
	}
	
	if(nypoint < 2) return;

	float yzsum =0;
	float z1sum = 0;
	float y2sum = 0;
	float ysum = 0;
	
	for(int ip=0; ip<nypoint; ip++){
		#pragma HLS loop_tripcount min=1 max=28
//		#pragma HLS pipeline II=1
		#pragma HLS unroll factor=1
		yzsum += (final_ypoints[ip].x * final_ypoints[ip].z);
		ysum += final_ypoints[ip].x;
		y2sum += (final_ypoints[ip].x * final_ypoints[ip].x);
		z1sum += final_ypoints[ip].z;
		}
	
	float my_initial = (nypoint*yzsum - ysum*z1sum)*1./max(nypoint*y2sum - ysum*ysum,float(1.e-3));
	float cy_initial = (z1sum - my_initial*ysum)*1./nypoint;
	
	// cout<<"my_inital "<<my_initial<<" cy_initial "<<cy_initial<<endl;
	
	float ychi2 = 0;
	for(int ip=0; ip<nypoint; ip++){
		#pragma HLS loop_tripcount min=1 max=28
//		#pragma HLS pipeline II=1
		#pragma HLS unroll factor=1
		float ypred = (my_initial*final_ypoints[ip].x + cy_initial);
		float yerr = my_initial*strip_width*1./sqrt(12);
		float chi2_try = (final_ypoints[ip].z - ypred)*(final_ypoints[ip].z - ypred)*1./(yerr*yerr);
		
		ychi2 += chi2_try;
		}	
	
	cout<<"ychi2/NDF:"<<ychi2<<"/"<<nypoint-2<<endl;
*/
	
/*

	for(int ip=0; ip<nplanes; ip++){
	
		float tryy = (zval[ip] - cy_initial)*1./my_initial;
		float min_dis = 1000;
		int npoints_plane = 0;
	
		for(int fp=0; fp<final_ypoints.size(); fp++){
			
			if(fabs(zval[ip] - final_ypoints[fp].z) < 1.e-6){
				
				float ydiff =  fabs(tryy - final_ypoints[fp].x);
				
				if(ydiff < min_dis) { min_dis = ydiff; }
				npoints_plane++;
			}
		}
		
		if(npoints_plane>1){
			for(int fp=0; fp<final_ypoints.size(); fp++){
				if(fabs(zval[ip] - final_ypoints[fp].z) < 1.e-6){
					float ydiff =  fabs(tryy - final_ypoints[fp].x);
					if(ydiff > min_dis){
						final_ypoints.erase(final_ypoints.begin()+fp);
						}
					}
				}
			}
		
	}
	
  float my_final, cy_final;	
	
  for(int iloop=0; iloop<10; iloop++){
	
	for(int ip=0; ip<nplanes; ip++){
		
		float tryy = (zval[ip] - cy_initial)*1./my_initial;
		
		for(int fp=0; fp<final_ypoints.size(); fp++){
			
			if(fabs(zval[ip] - final_ypoints[fp].z) < 1.e-6){
				
				float ydiff =  fabs(tryy - final_ypoints[fp].x);
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
	
	my_final = (final_ypoints.size()*yzsum - ysum*z1sum)*1./max(final_ypoints.size()*y2sum - ysum*ysum,float(1.e-3));
	cy_final = (z1sum - my_final*ysum)*1./final_ypoints.size();
	
	// cout<<"my_final "<<my_final<<" cy_final "<<cy_final<<endl;
	
	if((fabs(my_final - my_initial) < 0.2*fabs(my_final+0.001)) && (fabs(cy_final - cy_initial) < 0.2*fabs(cy_final+0.001))) break;
	
	my_initial = my_final;
	cy_initial = cy_final;
	
	}
	
	ychi2 = 0;
	for(int ip=0; ip<final_ypoints.size(); ip++){
			
		float ypred = (my_final*final_ypoints[ip].x+cy_final);
		float yerr = my_final*strip_width*1./sqrt(12);
		float chi2_try = (final_ypoints[ip].z - ypred)*(final_ypoints[ip].z - ypred)*1./(yerr*yerr);
	
		ychi2 += chi2_try; 
		}	
	
	// cout<<"chi2/NDF:"<<ychi2<<"/"<<final_ypoints.size()-2<<endl;
	
	
	cout<<"Number of points "<<nypoint<<endl;
	cout<<"And points are: "<<endl;
	for(int ip=0; ip<nypoint; ip++){
		#pragma HLS loop_tripcount min=1 max=28
		cout<<"point:"<<ip+1<<" y:"<<final_ypoints[ip].x<<" z:"<<final_ypoints[ip].z<<endl;
	}
	*/
	// y-z plane ends

//	return;
	
}
