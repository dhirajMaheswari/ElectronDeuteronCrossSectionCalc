//header file for the V_18 potential

#ifndef AV_18_H
#define AV_18_H
#endif

#define INVHBARC	0.197328 //this is inverse (hbar*c) in GeV

#include <iostream>
#include <math.h>
 using namespace std;
//prototypes of the functions that calculates the wf and partial waves
float fd_AV18(float , int);
float U_V18(float); //"s" partial wave

float W_V18(float); //"d" partial wave
 
float C1[12], D1[12], BM1[12]; //arrays declared

//definitions of the functions start here

float fd_AV18(float x, int i) {
  float ans;
 if(i==1) {
	C1[0] = 0.706699E+00; C1[1] = -0.169743E+00; C1[2] = 0.112368E+01; C1[3] = -0.852995E+01;
	C1[4] = 0.195033E+02; C1[5] = -0.757831E+02; C1[6] = 0.283739E+03;
	C1[7] = -0.694734E+03; C1[8] = 0.885257E+03; C1[9] = -0.720739E+03;
	C1[10] = 0.412969E+03; C1[11] = -0.103336E+03;

//	 C1[12] = { 0.706699e+00, -0.169743e+00, 0.112368e+01, -0.852995e+01, 0.195033e+02,-0.757831e+02, 0.283739e+03,
//	   -0.694734e+03, 0.885257e+03, -0.720739e+03, 0.412969e+03, -0.103336e+03 };
	
	D1[0] = 0.176655E-01;D1[1] =  -0.124551E+00; D1[2] = -0.108815E+01;D1[3] =  0.384848E+01;
	D1[4] =  -0.852442E+01;D1[5] =  0.209435E+02; D1[6] =  -0.490728E+02;
	D1[7] = 0.577382E+02; D1[8] = -0.127114E+01; D1[9] = -0.628361E+02;
	D1[10] =  0.581016E+02; D1[11] = -0.177062E+02;

//	 D1[12] = {0.176655e-01, -0.124551e+00, -0.108815e+01, 0.384848e+01, -0.852442e+01, 0.209435e+02, -0.490728e+02,
//	  0.577382e+02, -0.127114e+01, -0.628361e+02, 0.581016e+02, -0.177062e+02 };

	BM1[0] = 0.2316; BM1[1] = 1.0; BM1[2] = 1.5; BM1[3] = 2.0; BM1[4] = 2.5;
	BM1[5] =  3.5;BM1[6] =  4.5; BM1[7] = 5.5; BM1[8] = 6.5; BM1[9] = 8.0; BM1[10] = 9.5;BM1[11] =  11.0;

//	BM1[12] = {0.2316, 1.0, 1.5, 2.0, 2.5, 3.5, 4.5, 5.5, 6.5, 8.0, 9.5, 11.0 };
        ans = 0.0;
	return ans;
	}
 else{
   ans = (powf(U_V18(x/INVHBARC),2.0) + powf(W_V18(x/INVHBARC), 2.0))/powf(INVHBARC, 3.0);
  return ans;

  } 
} //end of the function

//"s" partial waves
 
 float U_V18(float x) {
  float A, F, ans;
  int j;	
  A = 0.0; F = 1.0;	
   for(j = 0;j<=11;j++) {
	A += C1[j]/(powf(x,2.0) + powf(BM1[j],2.0));
 } //end of for
 ans = A*F/ sqrt(4.0*M_PI);
 return ans;
} //end of function

//"d" partial waves

 float W_V18(float x) {
  float A, F, ans;
  int j;
  A = 0.0; F = 1.0;
  for(j = 0;j<=11;j++) {
    A += D1[j]/(powf(x,2.0) + powf(BM1[j],2.0));
  }
  ans = A*F/sqrt(4.0*M_PI);
 return ans;
} //end of funciton
