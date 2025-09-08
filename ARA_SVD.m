function [ U, S, V, H_w] = ARA_SVD( N,y,H,m_d, m_T)

% Differentiate y-vector
Y=zeros(N,1);
for i=2:N-1
    Y(i)=0.5*(y(i+1)-y(i-1));
end
Y(1)=0.5*y(2);
Y(N)=0.5*(y(N)-y(N-1));

% Weight matrix
if length(m_d) <= 1
    M   = sqrt(diag([m_d, 1, 1, 1, m_T]));
    % X*y all components (rho,u,v,w,T) for each grid point
    W   = kron(sqrt(diag(Y)),M);

else
    for jj = 1:length(m_d)
        M{jj}  = sqrt(Y(jj))*sqrt(diag([m_d(jj), 1, 1, 1, m_T(jj)]));
    end
    W  = blkdiag(M{:});
end


%% METHOD 1
% Operator for SVD
H_w = W*H*diag(1./diag(W));

% Filter reference Lele 1992 JCP
% q_vars = 5;
% H_w_filt = Filter_Operator(H_w,N,q_vars);

% H_w = H_w(1:1000,1:1000);
% W   = W(1:1000,1:1000);
% SVD gains
% [U0,S,V0] = svd(H_w);
[U0,S,V0] = svds(H_w,3);
% [U0,S,V0] = svds(H_w_filt,3);
% SIGMA = 5;
% options.subsppacedimension = 350;
% [U0,S,V0,flag] = svds(H_w,3,SIGMA,options);

U = diag(1./diag(W))*U0;
V = diag(1./diag(W))*V0;


% Normal check operator
% Norm = conj(U0)'*U0; % Norm = Id
% Normal check non-ortoghonal
% Norm = U_norm = (U)'*W*W*U; % Norm = Id

% figure; hold on
% for n_var = 2
%     Q{n_var}     = U(n_var:5:end,1); % Take all columns for each perturbation row
%     plot(y,Q{n_var})
% end
% plot(y,Q{4})
% plot(y(1:200),Q{2})


%% METHOD 2 - W_direct
% % Operator for SVD
% A_SVD = W*H*inv(W);
% % SVD gains
% [U,S,V] = svd(A_SVD);

%% METHOD 3 - TG
% % Operator for SVD
% work  = W*H;
% A_SVD = work'*work;
% % SVD gains
% [U,S,V] = svd(A_SVD);
% 
% S = sqrt(diag(S)); 




end