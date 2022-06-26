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
	GET_DM_VIEW(THETA, 2);

	GET_DM_VIEW(BTHETA, 3);
	GET_DM_VIEW(GAMMA, 2);
	GET_DM_VIEW(DELTA, 3);
	GET_DM_VIEW(varrho, 3);
	GET_DM_VIEW(zeta, 3);
	GET_DM_VIEW(BGAMMA, 3);
	GET_DM_VIEW(chi, 3);
	GET_DM_VIEW(ALPHAS, 2);

	// The following objects will be changed
	GET_DM_VIEW(C, 2);
	GET_DM_VIEW(D, 5);
	GET_DM_VIEW(E, 5);
	GET_DM_VIEW(G, 5);
	GET_DM_VIEW(H, 5);
	GET_DM_VIEW(M1, 3);
	GET_DM_VIEW(M2, 3);
	GET_DM_VIEW(M3, 3);

	GET_DM_VIEW(GD_sum, 4);
	GET_DM_VIEW(GC_sum, 4);
	GET_DM_VIEW(T_tilde, 3);

	GET_DM_VIEW(F1, 3);
	GET_DM_VIEW(F2, 3);
	GET_DM_VIEW(F3, 3);
	GET_DM_VIEW(F4, 3);

	GET_DM_VIEW(Q1, 3);
	GET_DM_VIEW(Q2, 3);
	GET_DM_VIEW(Q3, 3);
	GET_DM_VIEW(Q4, 3);

	
	
#pragma omp parallel for 
	for (int i = 1; i <= N; ++i) {
		for (int j = 1; j <= J ; ++j) {
			C(i, j) = -THETA(j,1) * GAMMA(j, i);
		}
	}


#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int i = 1; i <= N; ++i) {
			for (int j = 1; j <= J; ++j) {
				for (int n = 1; n <= N; ++n) {
					for (int m = 1; m <= N; ++m) {
						for (int o = 1; o <= N; ++o) {
							F(n, j, i, m, o, t) = THETA(j,1) * ((1.0- GAMMA(j, i))* DELTA(i + (j - 1)*N, o + (j - 1)*N, t) - DELTA(n + (j - 1)*N, o + (j - 1)*N, t)) * pi(o + (j - 1)*N, m, t);
							D(n, j, i, m, t) = D(n, j, i, m, t) - THETA(j, 1) * ((1.0 - GAMMA(j, i))*DELTA(i + (j - 1)*N, o + (j - 1)*N, t) - DELTA(n + (j - 1)*N, o + (j - 1)*N, t)) * pi(o + (j - 1)*N, m, t) * GAMMA(j, m);
							E(n, j, i, m, t) = E(n, j, i, m, t) + ((1.0 - GAMMA(j, i))* DELTA(i + (j - 1)*N, o + (j - 1)*N, t) - DELTA(n + (j - 1)*N, o + (j - 1)*N, t))* pi(o + (j - 1)*N, m, t);
						}
					}
				}
			}
		}
	}



#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int i = 1; i <= N; ++i) {
			for (int j = 1; j <= J; ++j) {
				for (int o = 1; o <= N; ++o) {
					for (int l = 1; l <= N; ++l) {
						G(i, j, o, l, t) = BGAMMA(i + (j - 1)*N, o + (j - 1)*N, t) * (1.0 - GAMMA(j, o))* varrho(l + (j - 1)*N, o, t);
					}
				}
			}
		}
	}


#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int i = 1; i <= N; ++i) {
			for (int j = 1; j <= J; ++j) {
				for (int o = 1; o <= N; ++o) {
					for (int k = 1; k <= J; ++k) {
						H(i, j, o, k, t) = BGAMMA(i + (j - 1)*N, o + (j - 1)*N, t) * (ALPHAS(j, o))* zeta(o + (k - 1)*N, j, t);
					}
				}
			}
		}
	}

#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int i = 1; i <= N; ++i) {
			for (int j = 1; j <= J; ++j) {
				for (int n = 1; n <= N; ++n) {
					//M1_temp(i + (j - 1)*N, i + (j - 1)*N, t, n) = GAMMA(j, i)*chi(n + (j - 1)*N, i, t) * C(i, j);
					M1(i + (j - 1)*N, i + (j - 1)*N, t) = M1(i + (j - 1)*N, i + (j - 1)*N, t) +GAMMA(j, i)*chi(n + (j - 1)*N, i, t) * C(i, j);
				}
			}
		}
	}


#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int n = 1; n <= N; ++n) {
			for (int j = 1; j <= J; ++j) {
				for (int m = 1; m <= N; ++m) {
					for (int o = 1; o <= N; ++o) {
						for (int l = 1; l <= N; ++l) {
							GD_sum(n, j, m, t) = GD_sum(n, j, m, t) + G(n, j, o, l, t)*D(l, j, o, m, t);
						}
					}
				}
			}
		}
	}

#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int j = 1; j <= J; ++j) {
			for (int m = 1; m <= N; ++m) {
				for (int i = 1; i <= N; ++i) {
					for (int n = 1; n <= N; ++n) {
						//M2_temp(i + (j - 1)*N, m + (j - 1)*N, t, n) = GAMMA(j, i)*chi(n + (j - 1)*N, i, t) * (D(n, j, i, m, t) + GD_sum(n, j, m, t));
						M2(i + (j - 1)*N, m + (j - 1)*N, t)= M2(i + (j - 1)*N, m + (j - 1)*N, t)+ GAMMA(j, i)*chi(n + (j - 1)*N, i, t) * (D(n, j, i, m, t) + GD_sum(n, j, m, t));
					}
				}
			}
		}
	}


#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int n = 1; n <= N; ++n) {
			for (int j = 1; j <= J; ++j) {
				for (int o = 1; o <= N; ++o) {
					for (int l = 1; l <= N; ++l) {
//						GC_temp(n, j, o, t, l) = G(n, j, o, l, t)*C(o, j);
						GC_sum(n, j, o, t) = GC_sum(n, j, o, t) + G(n, j, o, l, t)*C(o, j);
					}
				}
			}
		}
	}



#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int i = 1; i <= N; ++i) {
			for (int j = 1; j <= J; ++j) {
				for (int o = 1; o <= N; ++o) {
					for (int n = 1; n <= N; ++n) {
					//	M3_temp(i + (j - 1)*N, o + (j - 1)*N, t, n) = GAMMA(j, i)*chi(n + (j - 1)*N, i, t) * GC_sum(n, j, o, t);
						M3(i + (j - 1)*N, o + (j - 1)*N, t)= M3(i + (j - 1)*N, o + (j - 1)*N, t)+ GAMMA(j, i)*chi(n + (j - 1)*N, i, t) * GC_sum(n, j, o, t);
					}
				}
			}
		}
	}

#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int i = 1; i <= N; ++i) {
			for (int j = 1; j <= J; ++j) {
				for (int o = 1; o <= N; ++o) {
					for (int k = 1; k <= J; ++k) {
						for (int n = 1; n <= N; ++n) {
							//T_tilde_temp(i + (j - 1)*N, o + (k - 1)*N, t, n) = GAMMA(j, i)*chi(n + (j - 1)*N, i, t)*H(n, j, o, k, t);
							T_tilde(i + (j - 1)*N, o + (k - 1)*N, t)= T_tilde(i + (j - 1)*N, o + (k - 1)*N, t)+ GAMMA(j, i)*chi(n + (j - 1)*N, i, t)*H(n, j, o, k, t);
						}
					}
				}
			}
		}
	}


#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int i = 1; i <= N; ++i) {
			for (int j = 1; j <= J; ++j) {
				for (int n = 1; n <= N; ++n) {
					//Q1_temp(i + (j - 1)*N, i + (j - 1)*N, t, n) = GAMMA(j, i)*chi(n + (j - 1)*N, i, t);
					Q1(i + (j - 1)*N, i + (j - 1)*N, t) = Q1(i + (j - 1)*N, i + (j - 1)*N, t)+ GAMMA(j, i)*chi(n + (j - 1)*N, i, t);
					F1(i + (j - 1)*N, n + (i - 1)*N + (j - 1)*N*N, t) = -THETA(j,1)*GAMMA(j, i)*chi(n + (j - 1)*N, i, t); 
				}
			}
		}
	}


#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int i = 1; i <= N; ++i) {
			for (int j = 1; j <= J; ++j) {
				for (int m = 1; m <= N; ++m) {
					for (int h = 1; h <= N; ++h) {
						for (int n = 1; n <= N; ++n) {
							//F2_temp(i + (j - 1)*N, h + (m - 1)*N + (j - 1)*N*N, t, n) = -GAMMA(j, i)*chi(n + (j - 1)*N, i, t)*F(n, j, i, m, h, t);
							F2(i + (j - 1)*N, h + (m - 1)*N + (j - 1)*N*N, t) = F2(i + (j - 1)*N, h + (m - 1)*N + (j - 1)*N*N, t) -GAMMA(j, i)*chi(n + (j - 1)*N, i, t)*F(n,j,i,m,h,t);
						}
					}
				}
			}
		}
	}



	F(n, j, i, m, o, t)

#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int i = 1; i <= N; ++i) {
			for (int j = 1; j <= J; ++j) {
				for (int m = 1; m <= N; ++m) {
					for (int n = 1; n <= N; ++n) {
						//Q2_temp(i + (j - 1)*N, m + (j - 1)*N, t, n) = GAMMA(j, i)*chi(n + (j - 1)*N, i, t)*E(n, j, i, m, t);
						Q2(i + (j - 1)*N, m + (j - 1)*N, t)= Q2(i + (j - 1)*N, m + (j - 1)*N, t)+ GAMMA(j, i)*chi(n + (j - 1)*N, i, t)*E(n, j, i, m, t);
					}
				}
			}
		}
	}



#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int i = 1; i <= N; ++i) {
			for (int j = 1; j <= J; ++j) {
				for (int o = 1; o <= N; ++o) {
					for (int n = 1; n <= N; ++n) {
						for (int l = 1; l <= N; ++l) {
							Q3(i + (j - 1)*N, o + (j - 1)*N, t) = Q3(i + (j - 1)*N, o + (j - 1)*N, t) + GAMMA(j, i)*chi(n + (j - 1)*N, i, t)*G(n, j, o, l, t);
							F3(i + (j - 1)*N, l + (o - 1)*N + (j - 1)*N*N, t) = F3(i + (j - 1)*N, l + (o - 1)*N + (j - 1)*N*N, t) -THETA(j, 1)*GAMMA(j, i)*chi(n + (j - 1)*N, i, t)*G(n, j, o, l, t);
						}
					}
				}
			}
		}
	}



#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int i = 1; i <= N; ++i) {
			for (int j = 1; j <= J; ++j) {
				for (int m = 1; m <= N; ++m) {
					for (int n = 1; n <= N; ++n) {
						for (int o = 1; o <= N; ++o) {
							for (int l = 1; l <= N; ++l) {
								Q4(i + (j - 1)*N, m + (j - 1)*N, t)= Q4(i + (j - 1)*N, m + (j - 1)*N, t)+ GAMMA(j, i)*chi(n + (j - 1)*N, i, t)*G(n, j, o, l, t)*E(l, j, o, m, t);
							}
						}
					}
				}
			}
		}
	}


#pragma omp parallel for 
	for (int t = 1; t <= TIME; ++t) {
		for (int i = 1; i <= N; ++i) {
			for (int j = 1; j <= J; ++j) {
				for (int m = 1; m <= N; ++m) {
					for (int h = 1; h <= N; ++h) {
						for (int n = 1; n <= N; ++n) {
							for (int o = 1; o <= N; ++o) {
								for (int l = 1; l <= N; ++l) {
									F4(i + (j - 1)*N, h + (m - 1)*N + (j - 1)*N*N, t) = F4(i + (j - 1)*N, h + (m - 1)*N + (j - 1)*N*N, t) - GAMMA(j, i)*chi(n + (j - 1)*N, i, t)*G(n, j, o, l, t)*F(l, j, o, m, h, t);
								}
							}
						}
					}
				}
			}
		}
	}




}
 

