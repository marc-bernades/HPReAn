%% Load results and calculations
clear all; clc; close all;
% Select fluid and initial conditions
Fluid.Substance         = 'N2';
[~, T_c, P_c, ~, rho_c, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Fluid.Substance);
% Perturbation variables
q_vars = 5;
% Channel half height
delta  = 1;

% File name
label_save = 'HP_Turbulent'; % IG

%% Laminar baseflow
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_IG';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP_Re_10000'; 
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_IG_Re_10000'; 

%% Turbulent baseflow
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP_Turbulent';
name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP_Turbulent_PR';

% LSM
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP_Turbulent_PR_LSM';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP_Turbulent_PR_VLSM';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP_Turbulent_PR_NWSM_bw';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP_Turbulent_PR_NWSM_tw';

% Unpack structure and load variables
Results = load(strcat('Results/', name_file_load, '.mat'));
struct_field_names = fieldnames(Results.Data);
for k = 1:length(struct_field_names)
      eval([struct_field_names{k}, ' = Results.Data.',struct_field_names{k}, ';']);
end



%% Sigma 1 map - wavenumber space
if length(c) == 1
% Select Br and Re
Re_target = 1000;
[~,jj] = min(abs(n_LST_Sweep-Re_target));
Br_target = 5.6*1E-6;
[~,ii] = min(abs(N_target-Br_target));
Br_target = N_target(ii);
c_target = 0.88;
[~,cc] = min(abs(c-c_target));
c_target = c(cc);

for bb = 1:length(beta)
    for aa = 1:length(alpha)

        Sigma_1(aa,bb) = S{aa,bb,ii,jj,cc}(1,1);
    end
end

% Compute g
[alpha_grid,beta_grid]     = meshgrid(alpha,beta);
g = log10(Sigma_1'.^2.*alpha_grid.*beta_grid)';

% Apply log10
% Sigma_1     = log10(Sigma_1);
% Min and max adjustments for figure
Sigma_1_min = min(min(Sigma_1));
Sigma_1_max = max(max(Sigma_1)); %5000;
n_steps = 5;
Threshold = Sigma_1_max;
Sigma_1(Sigma_1 > Threshold) = Threshold;
Sigma_1(Sigma_1 < Sigma_1_min) = 0;

% Define log or linear
b_gain_log        = [0 4]; % 0 linear gain
b_gain_log        = [1 6]; % Re = 10000
b_gain_log        = [0 5]; % Turbulent

b_Ticks           = [1 2 3 4 5];
b_axis_labels_log = [-4 0; -2 2]; % 0 linear labels, [-1 4] log labels

plot_SigmaWavenumberSpace(alpha, beta, Br_target, Re_target, c, Sigma_1, ...
    Sigma_1_min, Sigma_1_max, n_steps, label_save, b_gain_log, b_axis_labels_log, b_Ticks)

end

%% Plot maximum singular vale at each dimension
% Select setups alpha and beta sweep RG vs IG
File_sweep{1} =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP';
File_sweep{2} =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_IG';

% Select setups c sweep RG vs IG
File_sweep{3} =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP_c_sweep';
File_sweep{4} =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_IG_c_sweep';

plot_TripletsEvol(File_sweep);

%% Family sigma values at single triplet
aa = 1; bb = 1; ii = 1; jj = 1; cc = 1;

% Singular value family at optimum kx = 1E-4 and kz = 2 at c = ub
% HP
name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP_singular_optima';

Results = load(strcat('Results/', name_file_load, '.mat'));
struct_field_names = fieldnames(Results.Data);
for k = 1:length(struct_field_names)
      eval([struct_field_names{k}, ' = Results.Data.',struct_field_names{k}, ';']);
end


Sigma_j_HP = diag(S{aa,bb,ii,jj,cc});

% IG
name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_IG_singular_optima';

Results = load(strcat('Results/', name_file_load, '.mat'));
struct_field_names = fieldnames(Results.Data);
for k = 1:length(struct_field_names)
      eval([struct_field_names{k}, ' = Results.Data.',struct_field_names{k}, ';']);
end

Sigma_j_IG = diag(S{aa,bb,ii,jj,cc});


f = figure
hold on; grid on; box on;
plot(Sigma_j_HP(1:20),'Marker','o','MarkerSize',8,'LineWidth',2,'LineStyle','-','color',[0.8500 0.3250 0.0980]);
plot(Sigma_j_IG(1:20),'Marker','o','MarkerSize',8,'LineWidth',2,'LineStyle','--','color',[0 0.4470 0.7410]);

xlabel('${j}$','interpreter','latex')
ylabel('${\sigma_{j}}$','interpreter','latex')

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)

set(gca,'YScale', 'log')
ylim([1 1E5])
yticks([1 10 100 1E3 1E4 1E5]); 

legend('$HP$','$IG$','interpreter','latex')
legend('Location','northeast','box','off')


pbaspect([1 0.8 1])

exportgraphics(f,strcat('Figures/','Sigma_family_HP_vs_IG','.jpeg'),'Resolution',300)



%% Analyse component-wise gains of H operator
H_comp_norm = zeros(length(alpha),length(beta),q_vars,q_vars);
H_norm      = zeros(length(alpha),length(beta));

for aa = 1:length(alpha)
    for bb = 1:length(beta)
        for n_var = 1:q_vars
            H_comp_norm(aa,bb,n_var,1)     = norm(H{aa,bb,ii,jj,cc}(n_var:q_vars:end,1));
            H_comp_norm(aa,bb,n_var,2)     = norm(H{aa,bb,ii,jj,cc}(n_var:q_vars:end,2));
            H_comp_norm(aa,bb,n_var,3)     = norm(H{aa,bb,ii,jj,cc}(n_var:q_vars:end,3));
            H_comp_norm(aa,bb,n_var,4)     = norm(H{aa,bb,ii,jj,cc}(n_var:q_vars:end,4));
            H_comp_norm(aa,bb,n_var,5)     = norm(H{aa,bb,ii,jj,cc}(n_var:q_vars:end,5));
        end
        H_norm(aa,bb)     = norm(H{aa,bb,ii,jj,cc});
    end
end

[beta_plot,alpha_plot]     = meshgrid(beta,alpha);
f = figure;
% Chccse component wise (1,1) rho rho, (1,2) rho u, ...(4,5) w T
% H_norm_plot = H_comp_norm(:,:,2,2);
% Or overall H operator
H_norm_plot = H_norm;
[c,h]=contourf(beta_plot,alpha_plot,log10(H_norm_plot'),'HandleVisibility','off');
set(h, 'edgecolor','none');
set(gca,'XScale', 'log')
set(gca,'YScale', 'log')
colorbar
colormap('jet')
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','log_{10}(\vert \vert H \vert \vert_{2})','$'],'Interpreter','latex','fontsize',14);
pbaspect([1 1.5 1])
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
xlabel('${\beta}$','interpreter','latex')
ylabel('${\alpha}$','interpreter','latex')


%% Plot responses and forcing
% Select Reynolds and Brinkman
Re_target = 1000;
[~,jj] = min(abs(n_LST_Sweep-Re_target));
Br_target = 5.6*1E-6;
[~,ii] = min(abs(N_target-Br_target));
Br_target = N_target(ii);
% Select c
c_target = 0.88;
[~,cc] = min(abs(c-c_target));
c_target = c(cc);
% Select wavenumbers
alpha_target = 0.0;
[~,aa] = min(abs(alpha-alpha_target));
alpha_target = alpha(aa);
beta_target = 2.0;
[~,bb] = min(abs(beta-beta_target));
beta_target = beta(bb);

% Select index
idx = 1; % Sigma 1

plot_ResponseForcing(y,delta,q_vars,U,S,V,aa,bb,ii,jj,cc,idx,alpha, beta, c, label_save)




%% Patterns
% Select Reynolds and Brinkman
Re_target = 1000;
[~,jj] = min(abs(n_LST_Sweep-Re_target));
Br_target = 5.6*1E-6;
[~,ii] = min(abs(N_target-Br_target));
Br_target = N_target(ii);
% Select c
c_target = 0.88;
[~,cc] = min(abs(c-c_target));
c_target = c(cc);
% Select wavenumbers
alpha_target = 0.0;
[~,aa] = min(abs(alpha-alpha_target));
alpha_target = alpha(aa);
beta_target = 2.0;
[~,bb] = min(abs(beta-beta_target));
beta_target = beta(bb);

% Select index
idx = 1; % Sigma 1

plot_PatternSigma(y,delta,q_vars,U,S,V,aa,bb,ii,jj,cc,idx,alpha, beta, c, label_save)


%% Scale motions
aa = 1; bb = 1; ii = 1; jj = 1; cc = 1;


% LSM PhD
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP_Turbulent_PR_LSM';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP_Turbulent_PR_VLSM';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP_Turbulent_PR_NWSM_bw';
% name_file_load =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP_Turbulent_PR_NWSM_tw';

% LSM CTR
idx       = 2; % Sigma 1/2
wall_plot = 'tw';
level                 = 6;   
name_file_load{1} =  'HP_Turbulent_VLSM';
name_file_load{2} =  'HP_Turbulent_LSM';
name_file_load{3} =  'HP_Turbulent_NWSM_bw';
name_file_load{4} =  'HP_Turbulent_NWSM_bw_classical';
name_file_load{5} =  'HP_Turbulent_NWSM_tw';
name_file_load{6} =  'HP_Turbulent_NWSM_tw_c5';
name_file_load{7} =  'HP_Turbulent_NWSM_tw_c10';
name_file_load{8} =  'HP_Turbulent_NWSM_tw_c12';

label_save   = {'VLSM','LSM','NWSM_bw','NWSM_bw_c','NWSM_tw','NWSM_tw_c5','NWSM_tw_c10','NWSM_tw_c12'}; % 'LSM', 'VLSM', 'NWSM_tw';
gamma_target = [pi,-pi/2,0,0,0,0,0,0]; % LSM = -pi/2, VLSM = pi, NWSM = 0

% Unpack structure and load variables
% Results = load(strcat('Results/PhD/', name_file_load, '.mat'));
Results = load(strcat('Results/', name_file_load{level}, '.mat'));

struct_field_names = fieldnames(Results.Data);
for k = 1:length(struct_field_names)
      eval([struct_field_names{k}, ' = Results.Data.',struct_field_names{k}, ';']);
end

% Outer Units
plot_ScaleMotion(y,delta,q_vars,U,S,V,aa,bb,ii,jj,cc,idx,alpha, beta, c, gamma_target(level),BF{1}, Dy{1}, label_save{level},wall_plot)
% Wall Units
plot_ScaleMotion_WallUnits(y,delta,q_vars,U,S,V,aa,bb,ii,jj,cc,idx,alpha, beta, c, gamma_target(level),BF{1}, Dy{1}, label_save{level},wall_plot)


