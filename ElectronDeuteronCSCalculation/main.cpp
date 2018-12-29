/*
	This code calculates the inclusive case inelastic cross section for the reaction e + d --> e' + X 
*/



#include "deuteronfn.cpp"
#include <fstream>
#include <iomanip>
#include <math.h>
#include "gadap.cpp"
#include "bodek.cpp"

#define massProton		0.938272
#define massNeutron		0.939565
#define massDeuteron	1.875628
#define pi				acos(-1.0)

float* edex(float, float, float, int, int, int);
/* first param: incident electron energy; second param: scattered electron energy; third param: electron scattering angle
	fourth param: Deuteron wf flag; fifth param: Bjorken limit flag; sixth param: indicates initialization or calculation

*/


float* sigmaInelastic(float, float, float, int, int, int);
/* first param: incident electron energy; second param: scattered electron energy; third param: electron scattering angle
   fourth param: initialize/PWIA; fifth param: Deuteron wf flag; sixth param: Bjorken limit

*/

float* polar_int(float, float, float);
/*
	arguments are the x,y,z components of momentum
	returns a pointer to array with p, theta and phi values

*/

float spectral_un(float, float, float, int, int); 

/*
	first param: recoil nucleon momentum; second param: recoil nucleon polar angle vs q
	third param: recoil nucleon azimuth angle; fourth param: who interacting: 1->proton, -1->neutron
	fifth param: initialization(0)/PWIA calculations(1)
	returns a float

*/


float* StructureFunctionNucleon(int, float, float, float,float);
/*
	first param: nucleon type; second param: effective x(x tilde); third param: q^2;
	fourth param: alphaN; fifth param: transverse momentum of interacting nucleon

	returns a pointer to two dimensional array such that first element is F1, second element is F2
*/


extern "C" {
		void struct_bodek_n_(float*,float*,float*,float*,float*);
	}

extern "C" {
		void struct_bodek_p_(float*,float*,float*,float*,float*);
	}


float thr_min(float);
float thr_max(float);

float undintInelastic(float, float);

//some global variables defined
float Ek, theta_e, Ekprime, BjorkenLimitFlag, xBj;
int ktr, int_case;

using namespace std;

int main() {

	int ix, flag, wfFlag;
	float* edexInit;
	float* edexCalc;
	float x, nu, QSq, dummy0, dummy1;
	fstream f;
	flag = 0; wfFlag = 1;
	

	edexInit = edex(Ek, Ekprime, theta_e, wfFlag, BjorkenLimitFlag, flag); // initialize
	cout<<"Init success"<<endl;
 

	Ek = 6.0; theta_e = 25.0*pi/180.0; //initial energy and scattering angle of electron11
	float th_half = theta_e/2.0;


//	f.open("data.dat",ios::out);
	for(ix = 0; ix <= 16; ix++) {
		x = 0.2 + float(ix)*0.05;
		xBj = x;
		dummy0 = massProton*x*Ek;
		dummy1 = massProton*x + 2.0*Ek*powf(sin(th_half), 2.0);
		Ekprime = dummy0/dummy1;
		nu = Ek - Ekprime;
		QSq = 4.0*Ek*Ekprime*powf(sin(th_half),2.0);
		x = QSq/(2.0*massProton*nu);
		xBj = x;
		flag = 1;	//calculation
			BjorkenLimitFlag = 0;
//		cout<<x<<"\t"<<QSq<<"\t"<<nu<<endl;


		edexCalc = edex(Ek, Ekprime, theta_e, wfFlag, BjorkenLimitFlag, flag); // calculations
		cout<<xBj<<"\t"<<QSq<<"\t"<<nu<<"\t"<<Ekprime<<"\t"<<edexCalc[8]<<endl;
//		f<<xBj<<"\t"<<QSq<<"\t"<<nu<<"\t"<<edexCalc[8]<<"\t"<<edexCalc[2]<<"\t"<<edexCalc[5]<<endl; //edexCalc[2] = W1total; edexCalc[5] = W2total;
	} 
		
//	f.close();
	return 0;
	
} //main ends	



float* edex(float ei, float epr, float thetae, int iw, int lbj, int icn) {
	/*	ei: initial energy GeV; epr: scattered electron's energy GeV; thetae: scattered electron's angle in radians;
		iw: flag to determine Deuteron wf (1: paris, 2: AV18); lbj: 0-no Bjorken limit, 1-yes Bjorken limit
		icn: 0 to initialize, 1 to calculate
		
		returns an array with nine (9) elements such that:
		first element: W1p, sec element: W1n, third element: W1sum;
		fourth element: W2p, fifth element: W2n, sixth element: W2sum;
		seventh element: crossSectionProtonInteracting; eight element: crossSectionNeutronInteracting; 
		ninth element: crossSectionTotal;
	
	*/

//	static float ans[9]; //array to hold final results


	float* sig_inInit;
	float* sig_inCalc;

	if(icn == 0) { //initialize
    	sig_inInit = sigmaInelastic(ei, epr, thetae, icn, iw,lbj);
		return sig_inInit;
	}
	
	else if(icn == 1) { //calculate
//		cout<<"got here."<<endl;
		sig_inCalc = sigmaInelastic(ei, epr, thetae, icn, iw,lbj);
		return sig_inCalc;
	}
	
 

} //end of edex function


float* sigmaInelastic(float ei, float epr, float thetae, int icn, int iw, int lbj) {
	/*	ei: initial energy GeV; epr: scattered electron's energy GeV; thetae: scattered electron's angle in radians;
		iw: flag to determine Deuteron wf (1: paris, 2: AV18); lbj: 0-no Bjorken limit, 1-yes Bjorken limit
		icn: 0 to initialize, 1 to calculate
		returns an array with nine (9) elements such that:
		first element: W1p, sec element: W1n, third element: W1sum;
		fourth element: W2p, fifth element: W2n, sixth element: W2sum;
		seventh element: crossSectionProtonInteracting; eight element: crossSectionNeutronInteracting; 
		ninth element: crossSectionTotal;	[crossSection is dsigma/(dOmega dE') in nb/(GeV.sr)
	*/
	
	static float ans[9];

	float pr, thr, phir, q2, q0, qv, x,y;
	float gmott;
//	int int_case;
	float pr_min, pr_max, eps, sum0;
	float W1proton, W1neutron, W1total;
	float W2proton, W2neutron, W2total;
	float csProton, csNeutron, csTotal;
	float  sum00;
	float F1D, F2D;
	float thetahalf = thetae/2.0;

	if(icn == 0) { //initialize
		float s_Un = spectral_un(pr, thr, phir, ktr, icn);

//		cout<<"Printing from sigmaInelastic function, init part."<<endl;
	}
	
	else if(icn == 1) { //calculations
//		cout<<"Calculations started for x:\t";
		q2 = 4.0*ei*epr*powf(sin(thetahalf),2.0); //Q^2 = 4EE'sin^2(theta/2)
		q0 = ei - epr; 
		qv = sqrt(q2 + q0*q0);
		x = q2/(2.0*massProton*q0);
		y = q0/ei;
//		cout<<q2<<"\t"<<q0<<"\t From sigmaInelastic function."<<endl;
		for(ktr = -1; ktr<=1; ktr+=2) {
			gmott = powf(ALPHA,2.0)*powf(cos(thetahalf),2.0)/(4.0*ei*ei*powf(sin(thetahalf),4.0))*0.389385*1000.0*1000.0;
			 int_case = 1;
			pr_min = 0.0; pr_max = 2.0; eps = 0.001;

//			cout<<ktr<<"   "<<q2<<"\t"<<q0<<"\t"<<x<<"\t"<<gmott<<endl;
		
			sum0 = gadap2(pr_min, pr_max, &thr_min, &thr_max, &undintInelastic, eps);
//			cout<<"Integral check.."<<endl;
			F1D = sum0*2.0*PI;

			 int_case = 2;
			pr_min = 0.0; pr_max = 2.0; eps = 0.001;
			sum00 = gadap2(pr_min, pr_max, &thr_min, &thr_max, &undintInelastic, eps);
			F2D = sum00*2.0*PI;
			
			if(ktr == 1) {
				W1neutron = F1D/massProton;
				W2neutron = F2D/q0; 	
			}
			else if(ktr == -1) {
				W1proton = F1D/massProton;
				W2proton = F2D/q0;
			}

		} //end ktr loop
		
		W1total = W1proton + W1neutron;
		W2total = W2proton + W2neutron;
		


		csNeutron = gmott*(W2neutron + 2.0*powf(tan(thetahalf), 2.0)*W1neutron);
		csProton  = gmott*(W2proton + 2.0*powf(tan(thetahalf), 2.0)*W1proton);  
		csTotal = (csProton + csNeutron);
		
//		csTotal = csTotal*y*y/(2.0*x*epr);// put in dsigma/dQ^2dx_A form;
		
		ans[0] = W1proton; ans[1] = W1neutron; ans[2] = W1total;
		ans[3] = W2proton; ans[4] = W2neutron; ans[5] = W2total;
		ans[6] = csProton; ans[7] = csNeutron; ans[8] = csTotal;

		return ans;
		
	}

	else {
		cout<<"Incorrect value passed."<<endl;
		exit(-1);
		}


}  //end of sigmaInelastic function 


float undintInelastic(float pr, float thr) {
	// this function gets integrated
	
	float phir, e_r, pr_z, pr_t, pst, ps, es, psz, pt;
	float alS, sinDelta, cosDelta;
	float q2, q0, qv;
	float alq, alphaN;
	float result;
	
	result = 0.0;
	phir = 0.0;
	e_r = sqrt(massProton*massProton + pr*pr);
	pr_z = pr*cos(thr);
	pr_t = pr*sin(thr); 
	
	alS = (e_r - pr_z)/massDeuteron*2.0;
//	cout<<pr<<"   "<<e_r<<"   "<<pr_z<<"   "<<pr_t<<"   "<<alS<<" From undint.. function"<<endl;

	if((alS <= 0.0) || (alS >= 2.0)) {
		return result; //just make sure alSpectator is in the range. it out off range, return 0,
	}
	

	float thetahalf = theta_e/2.0;

	q2 = 4.0*Ek*Ekprime*powf(sin(thetahalf),2.0); //Q^2 = 4EE'sin^2(theta/2)
	q0 = Ek - Ekprime; 
	qv = sqrt(q2 + q0*q0);
	
	
  
	sinDelta = sqrt(q2)/qv; cosDelta = q0/qv;
//	cout<<"q2: "<<q2<<" q0: "<<q0<<"  "<<sinDelta<<" "<<cosDelta<<"  "<<powf(sinDelta,2.0) + powf(cosDelta,2.0)<<" \t From undintInelastic function"<<endl;		

	alq = (q0 - qv)/massProton;	
	alphaN = 2.0 - alS; //momentum fraction of struck nucleon

	ps = pr; es = e_r;
	psz = pr_z; pt = pr_t;	
	//light cone variables definitions here
	
	float pplus = massDeuteron - (massProton*massProton + pt*pt)/(massProton*alS);
	float qplus = -q2/(massProton*alq);
	float pdotq = 0.5*(pplus*massProton*alq + qplus*massProton*alphaN);

	if(pdotq <= 0.0) return result;
//	cout<<q0<<"\t"<<pdotq<<endl;

	float xtil = q2/(2.0*pdotq);

	float starm2 = powf(massDeuteron - es, 2.0) - ps*ps;
	float w2nucleon = -q2 + 2.0*pdotq + starm2; //final produced mass squared
	float q0_off = (w2nucleon - massProton*massProton + q2)/(2.0*massProton);
	 xtil = q2/(2.0*massProton*q0_off);


//	cout<<q2<<"   "<<q0<<"   "<<xtil<<"\t from undint.. function"<<endl;
	if(BjorkenLimitFlag == 1) xtil = xBj/alphaN;


	if((xtil <= 0.0) || (xtil >= 1.0)) {

	
//	cout<<"   xtil: "<<xtil<<endl;
//		cout<<"xtil out off range."<<endl;
//		cout<<"alpha_s: "<<alS<<"\t"<<"x_tilde: "<<xtil<<endl;
		return result;
     }	

	int ktm = -ktr;
	float* F1F2Eff; //F1F2Eff[0] = F1eff, F1F2Eff[1] = F2Eff
	float F1Eff, F2Eff, F1D_siN, F2D_siN;
	
	// calculations of effective structure functions of nucleons

//	cout<<"Printing from undintIneastic funtion.."<<endl;
//	cout<<alS<<"\t"<<e_r<<"\t"<<pdotq<<"\t"<<BjorkenLimitFlag<<"\t"<<xtil<<endl;

	F1F2Eff = StructureFunctionNucleon(ktm, xtil, q2, alphaN, pt);
	

	F1Eff = F1F2Eff[0]; F2Eff = F1F2Eff[1];
	
	
	if(BjorkenLimitFlag == 0) {		
		float temp0 = powf(1.0 + cosDelta, 2.0); 
		float temp1 = powf(alphaN + pdotq/q2*alq, 2.0);	
		float temp2 = powf(sinDelta, 2.0)/2.0*powf(pt/massProton, 2.0);
		
		F2D_siN = (temp0*temp1 + temp2)*(q0*massProton)/pdotq*F2Eff;
		F1D_siN = F1Eff + powf(pt, 2.0)/(2.0*pdotq)*F2Eff;
		

//		F1D_siN = F1Eff;
//		 F2D_siN =0.0;// F2Eff; //for checking purposes

	}
	else if(BjorkenLimitFlag == 1) {
		F1D_siN = F1Eff;
		F2D_siN = alphaN*F2Eff;

	}
  
	
//	cout<<"F1eff:  "<<F1Eff<<"\t"<<"F2eff:  "<<F2Eff<<endl;

	// calculation of spectral function

	int icon = 1;
	float spec, Fsf;
	spec = spectral_un(pr, thr, phir, ktr, icon);

//	spec = 1.0;	//for checking purposes

	if(int_case == 1) {
		Fsf = F1D_siN*spec*pr*pr*sin(thr)/e_r; // d^3P_r/e_r is lorentz invriant phase = p_r^2dpr*sin(thr)d(thr)d(phi)/e_r
							
	}
	else if(int_case == 2) {
		Fsf = F2D_siN*spec*pr*pr*sin(thr)/e_r;
	}

	result = Fsf;

//	cout<<xBj<<"\t int_case:  "<<int_case<<"\t"<<"ktm: "<<ktm<<"\tFsf:  "<<Fsf<<endl;
	return result;
 

} //end undintInelastic function



float* StructureFunctionNucleon(int itn, float X, float Q2, float alN, float ptrans) {
/*
	This function returns the unpolarized structure functions of the nucleons using Bodek paramterization.
 		first param: identifies proton(1), neutron(-1)
		second param: Bjorken x value; third param: Q^2 value
		fourth param: alpha_N; momentum fraction of interacting nucleon
		fifth param: transverse momentum of the nucleon.
		return value: a pointer to two dim. array such that first element of array is F1, second element is F2.	
*/

	static float ans[2];

	float* f1f2;
	float F2a, F1a;
	float qzero = Q2/(2.0*X*massProton);
	float W2 = massProton*massProton + 2.0*qzero*massProton - Q2;
	
	if(itn == -1) {
//		f1f2 = strfn_n(qzero, Q2, W2);
		struct_bodek_n_(&qzero,&Q2,&W2,&F2a,&F1a);

	}	
	else if(itn == 1) {
//		f1f2 = strfn_p(qzero, Q2, W2);
		struct_bodek_p_(&qzero,&Q2, &W2, &F2a, &F1a);
	}
//	ans[0] = f1f2[0]; ans[1] = f1f2[1]; //ans[0] has f1, ans[1] has f2

	ans[0] = F1a; ans[1] = F2a;
//	cout<<"Printing from StructureFunctionNucleon function.."<<endl;
//	cout<<"f1: "<<ans[0]<<"   f2:  "<<ans[1]<<endl;
	
	return ans;

} //end of StructureFunctionNucleon function	


float spectral_un(float pr, float thr, float phir, int ktra, int icon) {
	/*
		This function calculates the spectral function of deuteron.
	*/
	
	int kj, ktm, ksm, ksr, iw, ini;
	float pmv, thm, phim, wfre, wfim;
	float pr_z, pr_t, pr_x, pr_y, als, pk, ek, pkz;
	float* polarIntRetVal;
	float ans, fct0, spec0;

	COMPLEX wfInit, wfCalc;

	
	iw = 1; //Paris potential
	if(icon == 0) {
		ini = 1; 
		wfInit = DeuteronWf(kj, ktm, ksm, ksr, pmv,  thm, phim, iw, ini);
	}

	else if(icon == 1) { //here, calculate the spectral function under PWIA

		float e_r = sqrt(massProton*massProton + pr*pr);
		pmv  = pr; thm = PI - thr; phim = PI + phir;
		fct0 = e_r;

		if(e_r < massDeuteron) fct0 = e_r*massDeuteron/(2.0*(massDeuteron - e_r));

		pr_z = pr*cos(thr);
		pr_t = pr*sin(thr);
		pr_x = pr_t*cos(phir); pr_y = pr_t*sin(phir);
		
		als = (e_r - pr_z)/massProton;
		float temp0 = powf(massProton,2.0) + powf(pr_t, 2.0);
		float temp1 = als*(2.0 - als);

		pk = sqrt(temp0/temp1 - massProton*massProton);
		ek = sqrt(massProton*massProton + pk*pk); 
		pkz = ek*(1.0 - als);
		polarIntRetVal = polar_int(pr_x, pr_y, pkz);
		pmv = polarIntRetVal[0]; 
		thm = PI - polarIntRetVal[1]; phim = PI + polarIntRetVal[2];

//		cout<<"alphaS:  "<<als<<"\t"<<pr_x<<"\t"<<pr_y<<"\t"<<pr_z<<"\t"<<pmv<<"\t"<<thm<<"\t"<<phim<<endl;
		fct0 = ek/(2.0 - als);

		fct0 = 1.0/(2.0-als);  // this with light-cone deuteron wf gives numbers which are closest

		ktm = -ktr;
		spec0 = 0.0; ini = 0;

		for(kj = -1; kj<=1; kj++) { 
		  for(ksm = -1; ksm <= 1; ksm+= 2) {
//			for(ksr = -1; ksm<= 1; ksm+= 2) { // this typo took two days to be figured out.
			for(ksr = -1; ksr<= 1; ksr+= 2) {
				wfCalc = LC_DeuteronWf(kj, ktm, ksm, ksr, als, pr_x, pr_y,iw, 0);
//				wfCalc = DeuteronWf(kj, ktm, ksm, ksr, pmv, thm, phim,iw, 0);
				spec0+= powf(wfCalc.real, 2.0) + powf(wfCalc.imag, 2.0);

			} //ksr loop ends
		  } //ksm loop ends
		} //kj loop ends

	} //icon = 1 block ends here
	
	ans = fct0*spec0/3.0;

//	cout<<"Result of spectral_un function:  "<<ans<<endl;
	return ans;

} //end of spectral_un function



float* polar_int(float px, float py, float pz) {

	float p = sqrt(px*px + py*py + pz*pz);
	float th, phi;
	static float ans[3];

	th = phi = 0.0;

	if(p == 0.0) {
		ans[0] = 0.0; ans[1] = 0.0; ans[2] = 0.0;
		return ans;
	}

		float arg = pz/p;
		if(arg > 1.0) arg = 1.0;
		if(arg < -1.0) arg = -1.0;
		th = acos(arg);
		float snth = sin(th);
		phi = 0.0;
		if(py == 0.0 && px < 0.0) { //a
			phi = acos(-1.0);
			ans[0] = p; ans[1] = th; ans[2] = phi;
			return ans;
		} //a
		else if(py == 0.0 && px > 0.0) { //b
			phi = 0.0;
			ans[0] = p; ans[1] = th; ans[2] = phi;
			return ans;
		} //b		
	

	 if(snth != 0.0) { //c
		 arg = px/(p*snth);
		if(arg > 1.0) arg = 1.0;
		if(arg < -1.0) arg = -1.0;
		phi = acos(arg);
		if(py/snth < 0.0) {
			phi = 2.0*acos(-1.0) - phi;
			}
			ans[0] = p; ans[1] = th; ans[2] = phi;
			return ans;
		} //c		
	
} //end of polar_int fnction



float thr_min(float q) {
	return 0.0;
}

float thr_max(float q) {
	return acos(-1.0);
}


