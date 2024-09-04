#include "armaMex.hpp"

int K_medians_target_function(mat X, rowvec a, vec Z){
	
	int N = X.n_rows;
	int K = X.n_cols;
	double KMsum = 0;
	
	for(int i=0; i<N; i++){
		KMsum += Z(i) * norm(X.row(i)-a, 2);
	}
	
	return KMsum;
	
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	// Check the number of input arguments.
	if (nrhs != 5)
		mexErrMsgTxt("Incorrect number of input arguments.");
	if (nlhs != 2)
		mexErrMsgTxt("Incorrect number of output arguments.");
	
	// Check type of input.
	//if ( (mxGetClassID(prhs[0]) != mxDOUBLE_CLASS) || (mxGetClassID(prhs[1]) != mxDOUBLE_CLASS) || (mxGetClassID(prhs[2]) != mxINT32_CLASS) || (mxGetClassID(prhs[3]) != mxINT32_CLASS) )
	//mexErrMsgTxt("Incorrect input type.");
	//
	// variable clain
	//
	
	//// read from input
	mat X = armaGetPr(prhs[0]);
	int K = mxGetScalar(prhs[1]);
	double step_size = mxGetScalar(prhs[2]);
	int round_MAX = mxGetScalar(prhs[3]);
	int grad_MAX = mxGetScalar(prhs[4]);
	
	//// define others
	int i, k, count_fail;
	int N = X.n_rows;
	int p = X.n_cols;
	
	uword whichmin;
	rowvec rowcopy(p);
	rowvec Zrowcopy(K);
	rowvec colsums(K);
	int n;
	int fixed_start;
	
	mat Z(N, K);  // Z: membership matrix
		// initialize Z
		int h = (int)(N/K);
		Z.zeros();
		
		// random start
		for(n=0; n<N; n++){
			Zrowcopy = randu<rowvec>(K);
			double min_val = Zrowcopy.min(whichmin);
			Z.row(n) = zeros(1, K);
			Z(n, whichmin) = 1;
		}
		// if random start fails to produce non-degenerated columns then fixed start:
		fixed_start = 0;
		colsums = sum(Z, 0);
		if( min(colsums)/N < 0.05 ){
			fixed_start = 1;
		}
		if(fixed_start==1){
			Z.col(K-1) = ones(N, 1);
			for(k=0; k<K-1; k++){
				Z.rows(k*h, (k+1)*h-1) = zeros(h, K);
				Z.submat(k*h, k, (k+1)*h-1, k) = ones(h, 1);
			}
		}
	mat C(K, p);  // C: community centers
		// initialize C
		for(k=0; k<K; k++){
			C.row(k) = sum((Z.col(k)*ones(1, p)) % X, 0)/sum(Z.col(k));
		}
		
	mat Distances(N, K);
	
	// for center updates
	rowvec gradk(p);
	mat Xk(N, p);
	double step;
	
	//
	// end variable claim
	//
	
	
	//
	// iteration
	//
	
	for(int round_ind=0; round_ind<round_MAX; round_ind++){
		
		// 1. update centers
		for(k=0; k<K; k++){
			
			// calculate gradient
			for(int grad_ind=0; grad_ind<grad_MAX; grad_ind++){
				Xk = ones(N, 1)*C.row(k) - X;
				gradk.zeros();
				for(i=0; i<N; i++){
					if( Z(i, k)*norm(Xk.row(i), 2)>1e-5 ){
						gradk += Xk.row(i)/sqrt(norm(Xk.row(i), 2));
					}
				}
			}
			
			// use gradient to update centers
			step = step_size;	// working step size		
			count_fail = 0;
			while( K_medians_target_function(X, C.row(k)-step*gradk, Z.col(k)) > K_medians_target_function(X, C.row(k), Z.col(k)) && count_fail<10){
				step *= 0.5;
			}
			if(count_fail<10){
				C.row(k) -= step*gradk;
				step *= 1.2;
			}
			// else: do not update
			
		}
		
		
		// 2. update grouping labels
		for(k=0; k<K; k++){
			Distances.col(k) = sqrt(sum( square(X - ones(N, 1)*C.row(k)) , 1));
		}
		Z.zeros();
		for(i=0; i<N; i++){
			rowcopy = Distances.row(i);
			rowcopy.min(whichmin);
			Z(i, whichmin) = 1;
		}
		// end update grouping labels
		
	}
	
	plhs[0] = armaCreateMxMatrix(C.n_rows, C.n_cols);
	armaSetPr(plhs[0], C);
	plhs[1] = armaCreateMxMatrix(Z.n_rows, Z.n_cols);
	armaSetPr(plhs[1], Z);
	
	/*return List::create(Named("centers")=C, Named("clusters")=Z);*/
	
}
