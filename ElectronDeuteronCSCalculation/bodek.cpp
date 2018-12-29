/*
This is the bodek parametrization which will be used to parametrize the 
unpolarized structure functions F1 and F2 for proton and neuteron.

*/

# define m_p	0.9382727 //GeV
# define PI		acos(-1.0)
# define ALPHA	1.0/137.0



float Maximum(float, float);
float Bodek(float, float);
float gp_h(float, float);
float gp_n(float, float);
float* strfn_p(float, float, float); // returns the f1 and f2 for proton
float* strfn_n(float, float, float); // returns f1 and f2 for neutron


  float gp_h(float q0, float Q2) {
	float xx, gi, ww, t, wp, res;
	float pm = 0.938279;
	xx = Q2/(2.0*pm*q0); 
	gi = 2.0*pm*q0;
	ww = (gi + 1.642)/(Q2 + 0.376);
	t = 1.0 - 1.0/ww;
	wp = 0.256*powf(t, 3.0) + 2.178*powf(t, 4.0) + 0.898*powf(t, 5.0) - 6.716*powf(t, 6.0) + 3.756*powf(t, 7.0);
	res = wp*ww*Q2/gi;
	return res;
 }

  float gp_n(float q0, float Q2) {
	float xx, gi, ww, t, wn, res;
	float pm = 0.938279;	
	xx = Q2/(2.0*pm*q0); 
	gi = 2.0*pm*q0;
	ww = (gi + 1.642)/(Q2 + 0.376);
	t = 1.0 - 1.0/ww;
	wn = 0.064*powf(t, 3.0) + 0.225*powf(t, 4.0) + 4.106*powf(t, 5.0) - 7.079*powf(t, 6.0) + 3.055*powf(t, 7.0);
	res = wn*ww*Q2/gi;
	return res;
 }


 float* strfn_p(float q0, float Q2, float w2) {
/*
	since C++ doesn't allow to return an array and it does not advocate to return the address of a local variable to outside of the function, so you would have to define the local variable as static variable. 

the first element of the returned array is F1 and second is F2.


*/
	

	float f1, f2, r, w2h, wi, w1h, pm;
	static float ans[2];
	pm = 0.938279; r = 0.18; w2h = 0.0;
	
	if(w2 < powf(pm, 2.0)){
		ans[0] = ans[1] = 0.0;
	 return ans;
	}
	wi = sqrt(w2); 
	w2h = gp_h(q0,Q2)*Bodek(wi, Q2)/q0;
	w1h = (1.0 + powf(q0, 2.0)/Q2)/(1.0 + r)*w2h;
	f2 = q0*w2h;
	f1 = pm*w1h;
	ans[0] = f1; ans[1] = f2;
	return ans;

 }


float* strfn_n(float q0, float Q2, float w2) {
	float f1, f2, r, w2nt, wi, w1nt, pm;
	static float ans[2];
	pm = 0.938279; r = 0.18; w2nt = 0.0;
	
	if(w2 < powf(pm, 2.0)) {
		ans[0] = ans[1] = 0.0;
	 return ans;
	}
	wi = sqrt(w2); 
	w2nt = gp_n(q0,Q2)*Bodek(wi, Q2)/q0;
	w1nt = (1.0 + powf(q0, 2.0)/Q2)/(1.0 + r)*w2nt;
	f2 = q0*w2nt;
	f1 = pm*w1nt;
	ans[0] = f1; ans[1] = f2;
	return ans;

 }




  float Bodek (float WM, float QSQ) {
	int LSPIN[4] = {1, 2, 3, 2};
	float C[24] = {1.0741163, 0.75531124, 3.3506491, 1.7447015, 3.5102405, 1.040004, 1.2299128, 0.10625394,
		       0.48132786, 1.5101467, 0.081661975, 0.65587179, 1.7176216, 0.12551987, 0.7473379, 1.953819,
			0.19891522, -0.17498537, 0.0096701919, -0.035256748, 3.5185207, -0.59993696,
			4.7615828, 0.41167589};
	float PM = 0.938256;
	float PM2 = 1.876512;
	float PMSQ = 0.880324;
	int NRES = 4;
	int NBKG = 5;
	float B = 0.0;
	float result, omega, X, XPX, PIEMSQ, B1, EB1, B2, EB2, WSQ;
	float BBKG, BRES, res, ressum, ram, rma, rwd, qstarn, qstar0, term, term0, brwig;
	float gamres;
	int i,j,k, index;

	res = 0.0;
	if(WM <= 0.939) return 0.0;
	WSQ = powf(WM, 2.0);
	omega = 1.0 + (WSQ - PMSQ)/QSQ;
	X = 1.0/omega;
	XPX = C[21] + C[22]*powf(X - C[23],  2.0);
	PIEMSQ = powf(C[0] - PM, 2.0);
	
	B1 = 0.0;

	if(WM == C[0]) goto aa;

	B1 = Maximum(0.0, WM - C[0])/(WM - C[0])*C[1];

	aa:
	    EB1 = C[2]*(WM - C[0]);
		
	  if(EB1 > 25.0) goto a;
	B1 = B1*(1.0 - exp(-EB1));
	B2 = 0.0;
	if(WM == C[3]) goto bb;
	a:
	  B2 = Maximum(0.0, WM - C[3])/(WM - C[3])*(1.0 - C[1]);
	  
	 
	bb:
	EB2 = C[4]*(WSQ - powf(C[3], 2.0));		
		
	if(EB2 > 25.0) goto ee; 
	 
  	B2 = B2*(1.0 - exp(-EB2));
	ee:	
	BBKG = B1 + B2;
	BRES = C[1] + B2;

	ressum = 0.0;
 	for( i = 0; i<NRES; i++) {
	  index = i*3 + 1 + NBKG;
//	cout<<index<<"\t";
	  ram = C[index];
	  if(i == 0) ram = C[index] + C[17]*QSQ + C[18]*powf(QSQ, 2.0);
	  rma = C[index + 1];
	  if(i == 2) rma = rma*(1.0 + C[19]/(1.0 + C[20]*QSQ));
	  rwd = C[index + 2];
	 qstarn = sqrt(Maximum(0.0, powf(((WSQ+PMSQ-PIEMSQ)/(2.*WM)),2.0) - PMSQ));
	  qstar0 = sqrt(Maximum(0.0, powf(((rma*rma-PMSQ+PIEMSQ)/(2.*rma)),2.0)-PIEMSQ));
	if(qstar0 <= 1.0e-10)  goto d;
	
	
	term = 6.08974*qstarn;
	term0 = 6.08974*qstar0;
	j = 2*LSPIN[i];
	k = j + 1;
	gamres = rwd*powf(term/term0, k)*(1.0 + powf(term0, j))/(1.0 + powf(term, j));
	 gamres = gamres/2.0;
	brwig = gamres/(PI*(powf(WM - rma, 2.0) + powf(gamres, 2.0)));
	res = ram*brwig/PM2;
	goto e;
	d:
	  res = 0.0;	
	e:
	 ressum+= res;
	}	

	result = BBKG*(1.0 + (1.0 - BBKG)*XPX) + ressum*(1.0 - BRES);
	
	return result; 

} //end of function


  float Maximum(float num1, float num2) {
     
    return((num1 > num2)?num1:num2);

 }

