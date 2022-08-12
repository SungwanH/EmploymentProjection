#include "mex.h"
#include "MatlabMatrix.h"
#include "stdio.h"
#include "CInterp.h"
#include "nr3_opt.h"

#include <iostream>
#include <boost/array.hpp>

#include <boost/numeric/odeint.hpp>
#include <omp.h>

using namespace std;



void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{	
	GET_INT(N_THREAD);
	omp_set_num_threads(N_THREAD);

	// The following objects will be read but not changed
	GET_INT(N);
	GET_INT(J);
	GET_INT(TIME);
	GET_DM_VIEW(pi, 3);
	GET_DM_VIEW(varrho, 3);
	GET_DM_VIEW(GAMMA, 2);

	// The following objects will be changed
	GET_DM_VIEW(BTHETA, 3);
	GET_DM_VIEW(U, 3);

#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int i = 1; i <= N; ++i) {
			for (int n = 1; n <= N; ++n) {
				for (int j = 1; j <= J ; ++j) {
					BTHETA(n + (j - 1)*N, i + (j - 1)*N, t) = pi(n + (j - 1)*N, i, t)*(1.0 - GAMMA(j, i));
				}
			}
		}
	}


#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int i = 1; i <= N; ++i) {
			for (int j = 1; j <= J; ++j) {
				for (int o = 1; o <= N; ++o) {
					U(i + (j - 1)*N, o + (j - 1)*N, t) = (1.0 - GAMMA(j, i))*varrho(o + (j - 1)*N, i, t);
				}
			}
		}
	}




}
 

