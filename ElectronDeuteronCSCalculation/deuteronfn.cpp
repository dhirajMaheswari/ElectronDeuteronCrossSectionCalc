#include "Paris.h"
#include "AV_18.h"
#include <stdlib.h>

#define INVHBARC	0.197328 //this is inverse (hbar*c) in GeV

 struct COMPLEX {
   float real;
   float imag; };

//functions prototypes

float uu(float, int); //s waves, int arg is a flag to check for the potential type used
float wd(float, int); //d waves, same defn for int flag
COMPLEX DeuteronWf(int, int, int, int, float, float, float, int,int); //calculates deuteron wavefunction
/* in this function, the first parameter takes values of 1, 0 ,-1 : deuteron spin projection
  second parameter takes values of 1,-1: struck out nucleon isospin [1:proton;-1:neutron]	
  third parameter {1,-1}: spin projection of proton
  fourth parameter {1,-1} : spin projection of neutron
  fifth param: momentum in GeV/c sixth param: polar angle in radians
 seventh param: azimuthal angle in radians
 eighth param: {1,2,3,4}::{Paris, V18, cdbonn, V14} used to select the type of potential 
ninth param: {1,0}::{1->initializing, 0->calculating}
 the function returns a complex number
*/

COMPLEX LC_DeuteronWF(int, int, int, int, float, float, float, int, int); //calculates the light cone deuteron wavefunction
/* in this function, the first parameter takes values of 1, 0 ,-1 : deuteron spin projection
  second parameter takes values of 1,-1: struck out nucleon isospin [1:proton;-1:neutron]	
  third parameter {1,-1}: spin projection of proton
  fourth parameter {1,-1} : spin projection of neutron
  fifth param: alpha, mom fraction ; sixth param: px component of momentum
 seventh param: py compoent of momentum
 eighth param: {1,2,3,4}::{Paris, V18, cdbonn, V14} used to select the type of potential 
ninth param: {0,1,2}::{0->initializing, 1-> z axis along total momentum of deuteron, 2 -> z axis opposite...}
 the function returns a complex number
*/
 

float MomentumDistribution(float); //calculates the momentum distribution
//functions definitions

float uu(float p, int potflag) {
float ans;
 if(potflag == 1) { //paris potential
  ans = U(p/0.197328)/sqrt(powf(0.197328,3.0));
 return ans;
   }
else if(potflag == 2) { //av18 potential
    ans = U_V18(p/0.197328)/sqrt(powf(0.197328,3.0));	
  return ans;
 }
} //end of uu

float wd(float p, int potflag) {
 float ans; 
 if(potflag == 1) {
 ans  = W(p/0.197328)/sqrt(powf(0.197328,3.0));
 return ans;
 }
 else if(potflag == 2) {
 ans  = W_V18(p/0.197328)/sqrt(powf(0.197328,3.0));
 return ans;
 }
  
}//end of wd

COMPLEX DeuteronWf(int kj, int kt, int ksp, int ksn, float p, float theta, float phi, int iw, int ini) {
   
  float wf_re, wf_im;
  float xxx;  
  float fis; //fis==1 implies proton was struck, -1 implies neutron was struck 
 COMPLEX result;
   wf_re = wf_im = 0.0;
   if(iw == 1) { //initialize the deuteron wf with specified potentials, 1 = paris, 2 = av18
	xxx = fd_Paris(0.0,1);
	}
   else if(iw == 2) {
   xxx = fd_AV18(0.0, 1);
    }	
     fis = 1.0; //assume proton was struck
   if(kt == -1) fis = -1.0; //neutron was struck
	
   if(kj == 1) { //A:   when deuteron spin projection is +1
 
	 //consider the case when both proton and neutron had their spins up	
	if((ksp == 1) && (ksn == 1)) { //B
	 wf_re = uu(p, iw) + wd(p, iw)/sqrt(8.0)*(3.0*powf(cos(theta),2.0)-1.0);
      wf_im = 0.0;
	} //end of B
  //consider the case when proton is up and neutron is down
 else if((ksp == 1) && (ksn == -1)) {//C
	wf_re = wd(p, iw)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*cos(phi);   
      wf_im = wd(p, iw)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*sin(phi);   
 	}//end of C
 //consider the case when proton is down and neutron is up
 else if((ksp ==-1) && (ksn == 1)) { //D
	wf_re = wd(p, iw)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*cos(phi);   
      wf_im = wd(p, iw)/sqrt(8.0)*3.0*cos(theta)*sin(theta)*sin(phi);   
      }//end of D
 //consider the case when proton and neutron both are down
 else if(ksp ==-1 && ksn ==-1) { //E
	wf_re = wd(p, iw)/sqrt(8.0)*3.0*powf(sin(theta),2.0)*cos(2.0*phi);
      wf_im = wd(p, iw)/sqrt(8.0)*3.0*powf(sin(theta),2.0)*sin(2.0*phi);
      }//end of E
    } //A	 
else if(kj ==0) {//A1: deuteron spin projection is 0
	//consider the case when both protn and neutron are up
	if((ksp ==1) && (ksn ==1)) { //B1
	wf_re = wd(p, iw) * 3.0 /2.0 * cos(theta)*sin(theta)*cos(phi);
	wf_im = - wd(p, iw) * 3.0 /2.0 * cos(theta)*sin(theta)*sin(phi);
	} //end of B1
	//consider the case when proton is up neutron is down
	else if((ksp == 1) && (ksn ==-1)) { //C1
	wf_re = uu(p, iw)/sqrt(2.0) - wd(p, iw)/2.0 * (3.0*powf(cos(theta),2.0) - 1.0);
	wf_im = 0.0;
	} //end of C
	//consider the case when proton is down neutron is up
	else if((ksp == -1) && (ksn == 1) ) { //D1
	wf_re = uu(p, iw)/sqrt(2.0) - wd(p, iw)/2.0*(3.0*powf(cos(theta),2.0) - 1.0);
	wf_im = 0.0;
	} //end of D1
	//when proton and neutron both are down
	else if(ksp == -1 && ksn == -1) { //E1
	wf_re = - wd(p, iw) * 3.0 /2.0 * cos(theta)*sin(theta)*cos(phi);
	wf_im = - wd(p, iw) * 3.0 /2.0 * cos(theta)*sin(theta)*sin(phi);
	} //end of E1
     } //end of A1

   else if(kj == -1) { //A2: deuteron spin projection is -1
	//consider when both spins are up
	if((ksp == 1) && (ksn == 1)) { //B2
	wf_re = wd(p, iw)/sqrt(8.0)*3.0*powf(sin(theta),2.0)*cos(2.0*phi);
	wf_im = - wd(p, iw)/sqrt(8.0)*3.0*powf(sin(theta),2.0)*sin(2.0*phi);
	} //end of B2
 	//proton up neutron down
    else if( (ksp == 1) && (ksn ==-1)) { //C2
	wf_re = - wd(p, iw)/sqrt(8.0) * 3.0*cos(theta)*sin(theta)*cos(phi);
	wf_im = wd(p, iw)/sqrt(8.0) * 3.0*cos(theta)*sin(theta)*sin(phi);
	}//end of C2
	//proton down neutron up
    else if((ksp == -1) && (ksn == 1)) { //D2
	wf_re = - wd(p, iw)/sqrt(8.0) * 3.0*cos(theta)*sin(theta)*cos(phi);
	wf_im = wd(p, iw)/sqrt(8.0) * 3.0*cos(theta)*sin(theta)*sin(phi);
	}// end of D2
   	//proton down neutron down
   else if( ksp == -1 && ksn == -1) { //E2
	wf_re = uu(p, iw) + wd(p, iw)/sqrt(8.0) * (3.0*powf(cos(theta),2.0) - 1.0);
	wf_im = 0.0;
	}//end of E2 
 	
} //end of A2



result.real = fis*wf_re;
result.imag = fis*wf_im;
 return result;
} //end of DeuteronWf function

float MomentumDistribution(float p) {
 COMPLEX WF;
 float res;
  float pi = acos(-1.0);
   int ini = 0; int kt = 1;
   int iw;	
   float theta = pi/10.0; float phi = pi/7.0;
   float sum = 0.0;
   int kj, ksp, ksn;
     for(kj = -1; kj<=1; kj++) { //1
	for(ksp = -1; ksp<=1; ksp+=2) { //2
	   for(ksn = -1; ksn <=1; ksn+=2) { //3
		WF = DeuteronWf(kj, kt, ksp, ksn, p, theta, phi, iw, ini);
		sum += powf(WF.real, 2.0) + powf(WF.imag, 2.0);       
	}//3		
     }//2
   } //1

 res = sum / 3.0; //averaged
 return res;
} //end of function


COMPLEX LC_DeuteronWf(int kj, int kt, int ksp, int ksn, float alpha, float px, float py, int iw, int inilc) {

	COMPLEX wf, result, f, f1;
	float pt, plc, elc, plcz, arg;
	float lcwf_real, lcwf_imag;
	 int ini;
	float mp = 0.938272; float p = 0.0; float theta = 0.0; float phi = 0.0;
	wf.real = 0.0; wf.imag = 0.0; result.real = 0.0; result.imag = 0.0;
	float mpn;
	if(kt == 1){
	 mpn = mp;
	}
	else if(kt == -1) {
	 mpn = 0.93956; //neutron mass
	}

	if(inilc == 0){
	  ini = 1;
	f = DeuteronWf(kj, kt, ksp, ksn, p,  theta, phi, iw, ini);
	} //end if

	if(alpha < 0.0 || alpha > 2.0) {
	  cout<<"alpha should be betwen 0 and 2."<<endl;
	  f1.real = -1.0; f1.imag = 0.0;
	  exit(-2);
	} //end alpha check

	//calculation part
	pt = sqrt(px*px + py*py);
	plc = sqrt((mpn*mpn + pt*pt)/(alpha*(2.0 - alpha)) - mpn*mpn);
	elc = sqrt(mpn*mpn + plc*plc);
	
	if(inilc == 1) {
	 plcz = (alpha - 1.0)*elc;	
	}
	else if(inilc == 2){
	 plcz = (1.0 - alpha)*elc;
	}

	//define angles
	  theta = 0.0;
	if(plc > 0.0){
	  arg = plcz/plc;
	  if(arg > 1.0) arg = 1.0;
	  if(arg < -1.0) arg = -1.0;	
	  theta = acos(arg);
	}

	phi = 0.0;
	if(pt > 0.0) {
	  arg = px/pt;
	  if(arg > 1.0) arg = 1.0;
	  if(arg < -1.0) arg = -1.0;
	  phi = acos(arg);
	 if(py < 0.0) phi = 2.0*M_PI - phi;
	}

	//calculate the wf
	
	 ini = 0;
	wf = DeuteronWf(kj, kt, ksp, ksn, plc, theta, phi, iw, ini);
	lcwf_real = wf.real*sqrt(elc);
	lcwf_imag = wf.imag*sqrt(elc);
	
	result.real = lcwf_real; result.imag = lcwf_imag;
	return result;	
	
} //end of function



