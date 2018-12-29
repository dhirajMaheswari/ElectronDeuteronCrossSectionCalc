
#include <iostream>
#include <math.h>

using namespace std;

float* polar_int(float, float, float);

int main() {

    float* ans;
	ans = polar_int(0.2, 0.2, 0.0);
	cout<<ans[0]<<"\t"<<ans[1]<<"\t"<<ans[2]<<endl;
 	return 0;

}
 


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
	
	} //end of function

