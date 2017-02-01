#include <math.h>
#include <mex.h>
#define PI 3.141592654

// simulate from exponential distribution
double exprnd(double rate, double p){
	return(-log(p)/rate);
}

void simulate_model(double *y_ret, double *times, double b1, double b2, int N, int num_obs, double *u){
	double t_curr = 0;
	int obs_counter=1, i, u_counter=0;
    int S_t = N;
    int I_t = 0;
    
	while(1){
        // time of next event
		t_curr += exprnd((b1 + b2*I_t)*S_t, u[u_counter]);
		u_counter++;
        // check if we can exit
        while (t_curr > times[obs_counter-1] && obs_counter <= num_obs){
            y_ret[obs_counter-1] = I_t;
            obs_counter++;
        }
		if (t_curr > times[num_obs-1]){ // then finished
			return;
		}				
		// now update the states
		I_t++; S_t--;	
        if (S_t==0){
            // then epidemic is OVER, fill in the rest of observations
            for (i = obs_counter; i <= num_obs; i++){
                y_ret[i-1] = N;
            }
            return;
        }
	}
	
}



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	int N, num_obs;
	double *times, *ret, b1, b2, *u;

	times = mxGetPr(prhs[0]);          // times to observe process
	b1 = mxGetScalar(prhs[1]);         // model parameter 1
    b2 = mxGetScalar(prhs[2]);         // model parameter 2
	N = (int)mxGetScalar(prhs[3]);     // initial number of susceptibles
	u = mxGetPr(prhs[4]);              // random numbers required to perform simulation
	num_obs = (int)mxGetNumberOfElements(prhs[0]);    // number of observations to observed
	
	plhs[0] = mxCreateDoubleMatrix(num_obs,1,mxREAL);
	
	ret = mxGetPr(plhs[0]);
	
	simulate_model(ret, times, b1, b2, N, num_obs, u);
	
	
}



