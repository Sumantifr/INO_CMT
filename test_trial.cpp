#include <stdio.h>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <typeinfo>
#include <vector>
#include <string.h>

int hextobin(string hex, int *out){

if(hex=="0") { out[3] = 0; out[2] = 0; out[1] = 0; out[0] = 0; }
//cout<<"0000"<<endl;	
if(hex=="1") { out[3] = 0; out[2] = 0; out[1] = 0; out[0] = 1; }
//cout<<"0001"<<endl;	
if(hex=="2") { out[3] = 0; out[2] = 0; out[1] = 1; out[0] = 0; } 
//cout<<"0010"<<endl;	
if(hex=="3") { out[3] = 0; out[2] = 0; out[1] = 1; out[0] = 1; } 
//cout<<"0011"<<endl;	
if(hex=="4") { out[3] = 0; out[2] = 1; out[1] = 0; out[0] = 0; } 
//cout<<"0100"<<endl;	
if(hex=="5") { out[3] = 0; out[2] = 1; out[1] = 0; out[0] = 1; } 
//cout<<"0101"<<endl;	
if(hex=="6") { out[3] = 0; out[2] = 1; out[1] = 1; out[0] = 0; } 
//cout<<"0110"<<endl;	
if(hex=="7") { out[3] = 0; out[2] = 1; out[1] = 1; out[0] = 1; } 
//cout<<"0111"<<endl;	
if(hex=="8") { out[3] = 1; out[2] = 0; out[1] = 0; out[0] = 0; } 
//cout<<"1000"<<endl;	
if(hex=="9") { out[3] = 1; out[2] = 0; out[1] = 0; out[0] = 1; } 
//cout<<"1001"<<endl;	
if(hex=="a") { out[3] = 1; out[2] = 0; out[1] = 1; out[0] = 0; } 
//cout<<"1010"<<endl;	
if(hex=="b") { out[3] = 1; out[2] = 0; out[1] = 1; out[0] = 1; } 
//cout<<"1011"<<endl;	
if(hex=="c") { out[3] = 1; out[2] = 1; out[1] = 0; out[0] = 0; } 
//cout<<"1100"<<endl;	
if(hex=="d") { out[3] = 1; out[2] = 1; out[1] = 0; out[0] = 1; } 
//cout<<"1101"<<endl;	
if(hex=="e") { out[3] = 1; out[2] = 1; out[1] = 1; out[0] = 0; } 
//cout<<"1110"<<endl;	
if(hex=="f") { out[3] = 1; out[2] = 1; out[1] = 1; out[0] = 1; }
//cout<<"1111"<<endl;	
return 0;
	
}

//int main(void)
void test_trial()
{
//	int i, k, l;
	std::ifstream file;
	file.open("eve.bin",std::ios::binary|std::ios::in);					// Reading the binary file
	
	long long int ref_res;
	const int nbyte = 2;
	unsigned char EventHeader[nbyte];
	int cnt = 0;

    bool istart = false; bool iend = false;
    
    ofstream fp_out;
    fp_out.open("INO_Data.txt");
    int ientry = 0;

	while(!(file.eof ()))
//	while(cnt<30)
	{
		file.read((char*) EventHeader, nbyte);
		unsigned int val, val1, val2, val3, val4;
		//val = (EventHeader[0]) |(EventHeader[1]<<8);////////(Reading the 4 digit block/16 bit hex word/2 layer info)
		val = (EventHeader[0]) | (EventHeader[1]<<8);
		val1 = (val)&0xf;
		val2 = (val>>4)&0xf;
		val3 = (val>>8)&0xf;
		val4 = (val>>12)&0xf;
		
		stringstream ss[4]; 
		ss[0] << hex << val1; 
		ss[1] << hex << val2; 
		ss[2] << hex << val3; 
		ss[3] << hex << val4; 
		
		std::cout<<std::hex<<val<<" -- "<<val1<<" -- "<<val2<<" -- "<<val3<<" -- "<<val4<<std::endl;
		
		int out1[4], out2[4], out3[4], out4[4];
		if(ss[0].str()=="a" && ss[1].str()=="a" && ss[2].str()=="a" && ss[3].str()=="a" ) { 
			istart = true; 
			ientry++;
			fp_out<<"Entry "<<ientry<<endl;
			continue;
			}
		if(ss[0].str()=="6" && ss[1].str()=="d" && ss[2].str()=="6" && ss[3].str()=="d" ) { istart = false; }
		if(istart){
			cout<<ss[0].str()<<" "<<ss[1].str()<<" "<<ss[2].str()<<" "<<ss[3].str()<<std::endl;
			}
		
		hextobin(ss[0].str(),out1);
		hextobin(ss[1].str(),out2);
		hextobin(ss[2].str(),out3);
		hextobin(ss[3].str(),out4);
//		cout<<out[3]<<" "<<out[2]<<" "<<out[1]<<" "<<out[0]<<endl;
		
		if(istart){
			
			fp_out<<out1[3]<<"\t"<<out1[2]<<"\t"<<out1[1]<<"\t"<<out1[0]<<"\t";
			fp_out<<out2[3]<<"\t"<<out2[2]<<"\t"<<out2[1]<<"\t"<<out2[0]<<endl;
			
			fp_out<<out3[3]<<"\t"<<out3[2]<<"\t"<<out3[1]<<"\t"<<out3[0]<<"\t";
			fp_out<<out4[3]<<"\t"<<out4[2]<<"\t"<<out4[1]<<"\t"<<out4[0]<<endl;
		}
		
		
		cnt++;
	}

fp_out.close();

}
