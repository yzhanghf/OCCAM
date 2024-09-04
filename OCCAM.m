function Z_hat = OCCAM(A, K, thresholding)
	% ASSUMING:
	% The network should be generated from assortative communities
	% INPUT:
	% A: adjacency matrix
	% K: number of communities
	% thresholding = 0,1: 1 means thresholding to discrete community assignments
	% OUTPUT:
	% Z_hat: estimated continuous community assignments
	% COPYRIGHT:
	% OCCAM: Overlapping Continuous Community Assignment Model
	% Reference: Detecting Overlapping Communities in Networks Using Spectral Methods
	% Yuan Zhang, Elizaveta Levina and Ji Zhu
	% Contact: Yuan Zhang yzhanghf@gmail.com
	
	N = size(A, 1);
	[U_A, S_A] = eig(A);
	positive = find(diag(S_A)>0);
	K1 = min(length(positive), K);
	pos = positive((length(positive)-K1+1):length(positive));
	U_A = U_A(:, pos);
	S_A = S_A(pos, pos);
	X_hat = U_A * sqrt(S_A);
		  
	X_tilde = X_hat;
	density_A = sum(sum(A))/(N*(N-1));
	tau = (density_A/K)^(1/5)*K^(3/2)/(N^(3/10)) / 10;
	for i=1:N
	X_tilde(i, :) = X_tilde (i, :) / (norm(X_tilde(i, :)) + tau);
	end
		  
	[centers, clusters] = K_medians_K_neq_p_MEX(X_tilde, K, 1e-2, 50, 50);
	optimal_centers = centers;
	optimal_clusters = clusters;
	optimal_loss =  sum( sqrt( sum( (optimal_clusters * optimal_centers - X_tilde).^2, 2 ) ) );
	for exp_ind = 1:9
		[centers, clusters] = K_medians_K_neq_p_MEX(X_tilde, K, 1e-2, 50, 50);
		loss = sum( sqrt( sum( (clusters * centers - X_tilde).^2, 2 ) ) );
		if(loss < optimal_loss)
			optimal_centers = centers;
			optimal_clusters = clusters;
			optimal_loss = loss;
		end
	end
	centers = optimal_centers;
	clusters = optimal_clusters;

	Z_tilde = X_tilde * transpose(centers) * pinv( centers*transpose(centers) );
	
	Z_hat = Z_tilde;
	Z_hat(Z_hat<0) = 0;
	for i=1:N
		Z_hat(i, :) = Z_hat (i, :) / (norm(Z_hat(i, :), 2) + tau/100);
		%%%% THRESHOLDING
		if(thresholding == 1)
			for k=1:K
				[z_max, z_which_max] = max(Z_hat(i, :));
				if(Z_hat(i, k)>threshold)
					Z_hat(i, k) = 1;
				else
					Z_hat(i, k) = 0;
				end
				% fail-safe measures
				if(sum(Z_hat(i, :))<0.001)
					Z_hat(i, z_which_max) = 1;
				end
				if(sum(Z_hat(i, :))<0.001)
					Z_hat(i, :) = 1;
				end
			end
		end
		%%%% END THRESHOLDING
	end
	
	
end % end function OCCAM