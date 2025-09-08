function H_filt = Filter_Operator(H,N,q_vars)

H_filt = zeros(size(H));

% Construct F operator
F = Filter_Matrix(N);

% Re-cast F discretization matrix
F_tot = zeros(N*q_vars,N*q_vars);
for ii = 1:length(F(:,1))
    vec_tot = [];
    for jj = 1:length(F(1,:))
        % Build the D ii-line layout given current q_vars
        vec_temp = zeros(1,q_vars);
        vec_temp(1) = F(ii,jj);
        vec_tot = [vec_tot vec_temp]; 
    end
    % Allocate pattern for all q_vars for D ii-line
    for n_var = 1:q_vars
        F_tot(ii*q_vars - (q_vars - 1) + n_var - 1,:) = circshift(vec_tot,n_var-1);
    end
end



H_filt = F_tot*H;
H_filt = F_tot*H*F_tot';





end