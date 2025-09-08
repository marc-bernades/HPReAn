% %% Linear Stability Theory Solver %%
% clear all; close all; clc

% Solver and model
bSolver           = 'Real';
HP_model          = 'HighPressure';  % CoolProp, RefProp, Constant, Power, LowPressure, HighPressure

% Fluid properties (ideal-gas)
Fluid.R_specific  = 188.9;
Fluid.gamma       = 1.289;
Fluid.mu_0        = 1;
Fluid.kappa_0     = 1;

%% Initial conditions
Fluid.Substance         = 'N2';
[~, T_c, P_c, ~, rho_c, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Fluid.Substance);
T_bw  = 0.75*T_c;
T_tw  = 1.5*T_c;
P_b   = 2.0*P_c;
% T_bw  = 310;
% T_tw  = 310;
% P_b   = 8E6;

% Grid points (y)
N = 200;
% N = 400; % For coherent structures

% Reference legnth scale
delta = 1; % Half channel-height
L_y   = 2*delta;
y_0   = 0; % Bottom y-position

%% Reference scaling
bScaling = 'bulk'; % wall, bulk

%% Grid
q_vars = 5;
% Chevysheb collocation Matrix
[y,D] = Chevysheb_collocation(N,delta); %dy = D*y = CentralDerivative_d1_2ndOrder_1D(y,y)

%% Baseflow input
Baseflow = 'Turbulent'; % Laminar (computed), Turbulent (DNS)

%% Baseflow
if strcmp(Baseflow,'Laminar')
    bTarget       = 'PrEc';
    N_target      = 5.6*1E-6;
    n_LST_Sweep   = 1000;

    BF        = cell(1,length(N_target));
    Dy        = cell(1,length(N_target));

    for ii = 1:length(N_target)

        disp("Computing Baseflow " + bTarget + " = "  + num2str(N_target(ii)) + ", Sweep " + num2str(ii) + " out of " + num2str(length(N_target)))
        [BF{ii}, Dy{ii}] = Baseflow_LST_adimensional(y,D,N,delta,L_y,T_bw,T_tw,P_b,bSolver,HP_model, Fluid.Substance, Fluid,bTarget,bScaling,N_target(ii));
        disp(" Resulted Reynolds " + bScaling + " = " + num2str(BF{ii}.Re) + newline)
    end

else

    % Turbulent mean flow (DNS)
    [BF, Dy] = Baseflow_Turbulent_mean_DNS(y,D,delta,bScaling); % Load DNS case
    % Same mimic targets
    bTarget       = 'PrEc';
    N_target      = BF{1}.Pr_b*BF{1}.Ec_b;
    n_LST_Sweep   = BF{1}.Re_b;
    

end

%% Adaptive resolvent analysis (ARA)
img   = sqrt(-1);

% Wavenumbers
% alpha = 1E-4; %logspace(-4,0,15); %0:0.2:2.0; % Streamwise wavenumber
% beta  = 2; %sort([logspace(-2,2,15) 2],'ascend'); %0:0.5:5.0; % Spanwise wavenumber
% c     = BF{1}.u_b./BF{1}.norm.u; % Phase speed
 
% Phase speed sweep
% alpha = 1E-4;
% beta  = 2;
% c     = BF{1}.u_b./BF{1}.norm.u; %logspace(-4,2,15); % Phase speed

% Coherent structures
% Plane of the velocity
% [~,u_y] = min(abs(BF{1}.y - y_target)); % y_target = 1, 1.9
% BF{1}.u_star(u_y);
% [~,u_y] = min(abs(BF{1}.y - 1.9));
% c_plus_pb = BF{1}.u_plus_tw(length(BF{1}.y) - u_y + 1);
% y_plus_pb = BF{1}.y_plus_tw(length(BF{1}.y) - u_y + 1);

% LSM
% lambda_x = 1; alpha = 2*pi/lambda_x;
% lambda_z = 0.5; beta  = 2*pi/lambda_z;
% c        = 0.95*BF{1}.u_b./BF{1}.norm.u; % Velocity center y = 1
% lambda_x_dim  = lambda_x*BF{1}.delta_h;
% lambda_x_plus = lambda_x_dim*(BF{1}.u_tau_bw/(BF{1}.mu_bw/BF{1}.rho_bw));
% lambda_z_dim  = lambda_z*BF{1}.delta_h;
% lambda_z_plus = lambda_z_dim*(BF{1}.u_tau_bw/(BF{1}.mu_bw/BF{1}.rho_bw));
% c_plus        = c*BF{1}.norm.u/BF{1}.u_tau_bw;
% [~,u_y]  = min(abs(BF{1}.u_plus_bw(1:100) - c_plus));
% y_plus   = BF{1}.y_plus_bw(u_y);

% VLSM
% lambda_x = 10; alpha = 2*pi/lambda_x;
% lambda_z = 1; beta  = 2*pi/lambda_z;
% c        = 1.1*BF{1}.u_b./BF{1}.norm.u; % Velocity of y_pb = 1.9
% lambda_x_dim  = lambda_x*BF{1}.delta_h;
% lambda_x_plus = lambda_x_dim*(BF{1}.u_tau_bw/(BF{1}.mu_bw/BF{1}.rho_bw));
% lambda_z_dim  = lambda_z*BF{1}.delta_h;
% lambda_z_plus = lambda_z_dim*(BF{1}.u_tau_bw/(BF{1}.mu_bw/BF{1}.rho_bw));
% c_plus        = c*BF{1}.norm.u/BF{1}.u_tau_bw;
% [~,u_y]  = min(abs(BF{1}.u_plus_bw(1:100) - c_plus));
% y_plus   = BF{1}.y_plus_bw(u_y);


% NWLSM BW
lambda_x_plus = 200;
lambda_x = lambda_x_plus/(BF{1}.u_tau_bw/(BF{1}.mu_bw/BF{1}.rho_bw));
lambda_x = lambda_x./BF{1}.delta_h; %Normalize DNS length
alpha    = 2*pi/lambda_x;
lambda_z_plus = 50;
lambda_z = lambda_z_plus/(BF{1}.u_tau_bw/(BF{1}.mu_bw/BF{1}.rho_bw));
lambda_z = lambda_z./BF{1}.delta_h; %Normalize DNS length
beta     = 2*pi/lambda_z;
c_plus   = 0.5;
c        = c_plus*BF{1}.u_tau_bw;
c        = c./BF{1}.norm.u;
[~,u_y]  = min(abs(BF{1}.u_plus_bw(1:200) - c_plus));
y_plus   = BF{1}.y_plus_bw(u_y);



% NWLSM TW
% lambda_x_plus = 450;
% lambda_x = lambda_x_plus/(BF{1}.u_tau_tw/(BF{1}.mu_tw/BF{1}.rho_tw));
% lambda_x = lambda_x./BF{1}.delta_h; %Normalize DNS length
% alpha    = 2*pi/lambda_x;
% lambda_z_plus = 100;
% lambda_z = lambda_z_plus/(BF{1}.u_tau_tw/(BF{1}.mu_tw/BF{1}.rho_tw));
% lambda_z = lambda_z./BF{1}.delta_h; %Normalize DNS length
% beta     = 2*pi/lambda_z;
% c_plus   = 10;
% c        = c_plus*BF{1}.u_tau_tw;
% c        = c./BF{1}.norm.u;
% [~,u_y]  = min(abs(BF{1}.u_plus_tw(1:100) - c_plus));
% y_plus   = BF{1}.y_plus_tw(u_y);
% y_norm   = 2 - BF{1}.y(u_y+1);

% Discounted resolvent - Spectrum check
b_eigenproblem = 0;
if b_eigenproblem == 1
    alpha = logspace(-4,0,15); %0:0.2:2.0; % Streamwise wavenumber
    beta  = sort([logspace(-2,2,15) 2],'ascend'); %0:0.5:5.0; % Spanwise wavenumber
    c     = BF{ii}.u_b./BF{1}.norm.u; % Phase speed
end

% Initialize
A            = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
B            = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep));
H            = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep),length(c));
H_w          = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep),length(c));
U            = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep),length(c));
V            = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep),length(c));
S            = cell(length(alpha),length(beta),length(N_target),length(n_LST_Sweep),length(c));
Sigma        = zeros(3,length(c));

for bb = 1:length(beta)
    for aa = 1:length(alpha)
        for ii = 1:length(N_target)
            for jj = 1:length(n_LST_Sweep)

                k = sqrt(alpha(aa).^2 + beta(bb).^2);

                disp("Computing LST " + "beta" + " = "  + num2str(beta(bb))  + " and " + "alpha" + " = "  + num2str(alpha(aa)) + " and " + bTarget + " = "  + num2str(N_target(ii)) + " Re = " + num2str(n_LST_Sweep(jj)) + ", Sweep " + num2str((bb-1)*length(alpha)*length(N_target)*length(n_LST_Sweep) + (aa-1)*length(N_target)*length(n_LST_Sweep) + (ii-1)*length(n_LST_Sweep) + jj) + " out of " + num2str(length(beta)*length(alpha)*length(N_target)*length(n_LST_Sweep)))

                %% Linear Operator (set omega = 0) as not needed here
                [A{aa,bb,ii,jj},B{aa,bb,ii,jj}] = LST(delta,N,BF{ii},Dy{ii},0,beta(bb),alpha(aa),bSolver,HP_model, Fluid.Substance, Fluid, bScaling, n_LST_Sweep(jj));

                if b_eigenproblem == 1
                % Eigenspectrum
%                 [Vec_eigen{aa,bb,ii,jj}, D_Val_eigen{aa,bb,ii,jj}] = eig(A{aa,bb,ii,jj},B{aa,bb,ii,jj},'qz');
%                 Val_eigen{aa,bb,ii,jj}                             = eig(A{aa,bb,ii,jj},B{aa,bb,ii,jj},'qz');

                N_modes    = 350; %length(y)*q_vars; % Modes to compute
                opts.v0    = rand(length(y)*q_vars,1);
                opts.tol   = 1e-15;
                opts.maxit = 1E10;
                sigma      = 0; % Shift value (e.g., around which you want to find eigenvalues)
                [Vec_eigen{aa,bb,ii,jj}, ~, ~] = eigs(A{aa,bb,ii,jj},B{aa,bb,ii,jj},N_modes,sigma,opts);
                Val_eigen{aa,bb,ii,jj}         = eigs(A{aa,bb,ii,jj},B{aa,bb,ii,jj},N_modes,sigma,opts);

                % Spectrum check
                for n_var = 1:q_vars
                    Val_eigen_var{n_var} = Val_eigen{aa,bb,ii,jj}(n_var:q_vars:end);
                    Vec_eigen_var{n_var} = Vec_eigen{aa,bb,ii,jj}(n_var:q_vars:end,:); % Take all columns for each perturbation row
                end
                realc = real(Val_eigen{aa,bb,ii,jj});
                imagc = imag(Val_eigen{aa,bb,ii,jj});
                Cond   = (abs(realc) <= 10) & (realc < 0.5) & imagc < 0.1;
                Val_eigen_target   = max(imagc(Cond));
                disp(['Most unstable mode = ', num2str(Val_eigen_target), ' total of = ', num2str(sum(imagc(Cond)>0))])

                % Plots
%                 x_lim = [-0.1 0.1];
%                 y_lim = [-0.5 0.1];
%                 plot_Spectrum_ARA(y, delta, N, alpha, N_target,n_LST_Sweep, Val_eigen{aa,bb,ii,jj}(Cond), 'test',x_lim,y_lim)
%                 plot_Perturbation(y, delta, N, alpha(aa), N_target(ii),  n_LST_Sweep(jj), Val_eigen{aa,bb,ii,jj}, Val_eigen_target, Vec_eigen_var, 'test')

                end


                for cc = 1:length(c)


                    omega = c(cc)*alpha(aa);
                   
                    disp("Computing resolvent " + num2str(cc) + " out of " + num2str(length(c)) + " at omega = " + num2str(omega))
                    disp(' ')
                    %% Resolvent Operator
%                     Id = eye(size(A{aa,bb,ii,jj}));
%                     Id(2:5,:) = 0; Id(end-3:end,:) = 0;
                    if b_eigenproblem == 1 && Val_eigen_target > 0
                        H{aa,bb,ii,jj,cc} = inv(A{aa,bb,ii,jj} + (- img.*omega + Val_eigen_target)*B{aa,bb,ii,jj}); %*Id;
                    else
                        H{aa,bb,ii,jj,cc} = inv(A{aa,bb,ii,jj} - img.*omega*B{aa,bb,ii,jj}); %*Id;
%                         lambda = eigs(H{aa,bb,ii,jj},N_modes,sigma,opts);
%                         Id = eye(size(A{aa,bb,ii,jj}));
%                         Id(2:5,:) = 0; Id(end-3:end,:) = 0;
%                         LAMBDA = img*A{aa,bb,ii,jj}*inv(B{aa,bb,ii,jj});
%                         lambda = diag(LAMBDA);
%                         plot_Spectrum_ARA(y, delta, N, alpha, N_target,n_LST_Sweep, lambda, 'test',x_lim,y_lim)

                    end

                    %% Resolvent analysis
                    % Weightings
                    m_d = 1;
                    m_T = 1;

                    [U{aa,bb,ii,jj,cc}, S{aa,bb,ii,jj,cc}, V{aa,bb,ii,jj,cc}, H_w{aa,bb,ii,jj,cc}] = ARA_SVD(length(y),y,H{aa,bb,ii,jj,cc},m_d, m_T);

                    N_modes    = 350; %length(y)*q_vars; % Modes to compute
                    opts.v0    = rand(length(y)*q_vars,1);
                    opts.tol   = 1e-15;
                    opts.maxit = 1E10;
                    sigma      = 0; % Shift value (e.g., around which you want to find eigenvalues)
                    Lambda     = eigs(A{aa,bb,ii,jj},B{aa,bb,ii,jj},N_modes,sigma,opts);
%                     Lambda     = eig(A{aa,bb,ii,jj},B{aa,bb,ii,jj});

                    plot_Spectrum_ARA(y, delta, N, alpha, N_target,n_LST_Sweep, Lambda, 'test',x_lim,y_lim)


                end
            end
        end
    end
end

%% Save results
name_label = 'HP_Turbulent_NWSM_bw_classical';


ARA_output(Fluid, T_bw, T_tw, T_c, P_b, P_c, N, y, D, delta, bScaling, bTarget,N_target, n_LST_Sweep, alpha, beta, c, BF,Dy,H, H_w, U, S, V, name_label);