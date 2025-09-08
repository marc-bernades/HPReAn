%% Load results and calculations
clear all; clc; close all;
% Select fluid and initial conditions
Fluid.Substance         = 'N2';
[~, T_c, P_c, ~, rho_c, ~,  ~,  ~,  ~,  ~,  ~, ~,  ~,  ~,  ~,  ~,  ~,  ~] = Substance_Library(Fluid.Substance);

%% File name
name_file_load = 'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_ARA_HP_Sweep';
% name_file_load = 'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_ARA_HP_Sweep_Ideal';


% Unpack structure and load variables
Results = load(strcat('Results/', name_file_load, '.mat'));
struct_field_names = fieldnames(Results.Data);
for k = 1:length(struct_field_names)
      eval([struct_field_names{k}, ' = Results.Data.',struct_field_names{k}, ';']);
end

% Check isothermal case
bIsothermal = contains(name_file_load,'Isothermal'); 
% Select Br to plot I- All, NI- only up to Br <= 0.1
N_plot = N_target;
% Perturbation variables
q_vars = 5;
% Channel half height
delta  = 1;

%% Map Transient Unstable points
GG_unstable = zeros(size(G_max));
for ii = 1:length(N_target)
    for aa = 1:length(alpha)
        for jj = 1:length(n_LST_Sweep)
            for bb = 1:length(beta)
                GG_gradient = gradient(GG{aa,bb,ii,jj});
                GG_smooth   = smooth(GG{aa,bb,ii,jj},10);
                Cond_max    = islocalmax(GG_smooth);
                Cond_unstable = max(GG_smooth(Cond_max)) < GG_smooth(end);
                if isempty(Cond_unstable) %No maxmimum
                    Cond_unstable = true;
                end
                if (max(islocalmax(GG{aa,bb,ii,jj})) == 0) || Cond_unstable
                    GG_unstable(aa,bb,ii,jj) = 1;
                end
            end
        end
    end
end

%% Transient growth
beta_target  = 2.0;
alpha_target = 0.0;
Re_target = 4500;
Br_target = 0.1;
[~,aa] = min(abs(alpha-alpha_target));
[~,bb] = min(abs(beta-beta_target));
[~,ii] = min(abs(N_target-Br_target));
[~,jj] = min(abs(n_LST_Sweep-Re_target));

% Sweep beta
figure
title(strcat('$','Re = ', num2str(n_LST_Sweep(jj)), ', Br = ',  num2str(N_target(ii)), ', \alpha = ',  num2str(alpha(aa)),'$'),'interpreter','latex')
for bb = 1:length(beta)
    plot(t_all{aa,bb,ii,jj}, GG{aa,bb,ii,jj}, 'DisplayName', strcat('$','\beta = ', num2str(beta(bb)), '$')); hold on
end

xlabel('${t}$','interpreter','latex')
ylabel('${G}$','interpreter','latex')

legend('interpreter','latex','Location','northeast','box','off')

set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

% Sweep Reynolds
figure
for jj = 1:length(n_LST_Sweep)
    plot(t_all{aa,bb,ii,jj}, GG{aa,bb,ii,jj}, 'DisplayName', strcat('$','Re = ', num2str(n_LST_Sweep(jj)), '$')); hold on
end

xlabel('${t}$','interpreter','latex')
ylabel('${G}$','interpreter','latex')

legend('interpreter','latex','Location','northeast','box','off')

set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

return

%% Modal counter-part
for n_var = 1:q_vars
    Val_eigen_var{n_var} = Val_eigen{aa,bb,ii,jj}(n_var:q_vars:end);
    Vec_eigen_var{n_var} = Vec_eigen{aa,bb,ii,jj}(n_var:q_vars:end,:); % Take all columns for each perturbation row
end

% Spectrum at given alpha and Re
plot_Spectrum(y, delta, N, alpha(aa), N_target(ii), n_LST_Sweep(jj), Val_eigen{aa,bb,ii,jj},  Vec_eigen_var, name_file_load)

% Select Br and Beta
if length(n_LST_Sweep) > 1
    beta_target = 0;
    Br_target   = 0.00001;
    [~,ii] = min(abs(N_target-Br_target));
    [~,bb] = min(abs(beta-beta_target));
    plot_StabilityDiagram(alpha, N_target, n_LST_Sweep, squeeze(GG_max(:,bb,:,:)), name_file_load, bIsothermal)
end


%% Perturbation at given mode
realc = real(Val_eigen{aa,bb,ii,jj});
imagc = imag(Val_eigen{aa,bb,ii,jj});
Cond   = (abs(realc) <= 10) & (realc < 0.5);
Val_eigen_target   = max(imagc(Cond));
plot_Perturbation(y, delta, N, alpha(aa), N_target(ii),  n_LST_Sweep(jj), Val_eigen{aa,bb,ii,jj}, Val_eigen_target, Vec_eigen_var, name_file_load)
plot_Perturbation_inset(y, delta, N, alpha(aa), N_target(ii),  n_LST_Sweep(jj), Val_eigen{aa,ii,jj}, Val_eigen_target, Vec_eigen_var, name_file_load)


%% Growth rate diagram
% Select Br and Re
Re_target = 1000;
[~,jj] = min(abs(n_LST_Sweep-Re_target));
% [~,ii] = min(abs(N_target-Br_target));

for ii = 1:length(N_target)
    if ~isempty(find(N_target(ii) == N_plot))
        Br_target = N_target(ii);
        GR = squeeze(GG_max(:,:,ii,jj));
        GR_min = 0;
        GR_max = 3000; % HP
%         GR_max = 400; % Ideal-gas
        Threshold = GR_max;
        GR(GR > Threshold) = Threshold;
        GR(GR < GR_min) = 0;
        n_steps = 5; %Plot colorbar
        % Set GR = 0 for alpha = 0, beta = 0
        % if alpha(aa) == 0 && beta(bb) == 0
        %     GR(1,1) = 0;
        % end
        plot_GrowthRate(alpha, beta, Br_target, Re_target, GR, name_file_load, GR_min, GR_max,n_steps)
    end
end

%% Optimal perturbation
beta_target  = 2;
alpha_target = 0;
Re_target = 1000;
Br_target = 0.50;
[~,aa] = min(abs(alpha-alpha_target));
[~,bb] = min(abs(beta-beta_target));
[~,ii] = min(abs(N_target-Br_target));
[~,jj] = min(abs(n_LST_Sweep-Re_target));

pos_in  = 4; % Normalize w
pos_out = 1; % Normalize u post_out = 2;

% Optimal perturbation and response
plot_OptimalPerturbation(BF{ii}.y, delta, q_vars, q_in{aa,bb,ii,jj}, q_out{aa,bb,ii,jj}, pos_in, pos_out,name_file_load, N_target(ii), alpha(aa), beta(bb))

%% Optimal perturbation pattern
beta_target  = 2;
alpha_target = 0;
Re_target = 1000;
Br_target = 5.6*1E-6;
[~,aa] = min(abs(alpha-alpha_target));
[~,bb] = min(abs(beta-beta_target));
[~,ii] = min(abs(N_target-Br_target));
[~,jj] = min(abs(n_LST_Sweep-Re_target));

pos_in  = 4; % Normalize w
pos_out = 2; % Normalize u

% Optimal perturbation pattern
for n_var = 1:q_vars
    Vec_eigen_var{n_var} = q_in{aa,bb,ii,jj}(n_var:q_vars:end); % Take all columns for each perturbation row
end

Norm = max(abs(Vec_eigen_var{pos_in}));
z      = linspace(0,1,200); % = Beta Z / (2pi)
% ph     = 0.45; % HP
ph     = -0.10; % IG

v      = Vec_eigen_var{3}/Norm; 
w      = Vec_eigen_var{4}/Norm;

% Optimal response pattern
for n_var = 1:q_vars
    Vec_eigen_var{n_var} = q_out{aa,bb,ii,jj}(n_var:q_vars:end); % Take all columns for each perturbation row
end

Norm = max(abs(Vec_eigen_var{pos_out}));

rho = Vec_eigen_var{1}/Norm; 
u   = Vec_eigen_var{2}/Norm; 
T   = Vec_eigen_var{5}/Norm; 

% Plot Spanwise (alpa = 0)
plot_OptimalPerturbationPattern(z,y,v,w,rho,u,T,ph,name_file_load,N_target(ii),alpha(aa),beta(bb))



%% Optimal perturbation multiple sets
% % Non-Isothermal flow cases (N4-5)
setup{1} = 'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_ARA_HP_Sweep';
setup{2} = 'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_ARA_HP_Sweep_Ideal';

setup_name = {'HP','IG'};
name_save  = 'HP_vs_IG';

% This should be based on current setup for input positions
beta_target  = 2;
alpha_target = 0;
Re_target = 1000;
% Br_target = 0.50; % Isothermal
Br_target = 5.6*1E-6; % Non-isothermal

pos_in  = 4; % Normalize w
% pos_out = 2; % Normalize u isothermal
pos_out = 1; % Normalize rho Non-isothermal

plot_OptimalPerturbation_MultipleCases(setup, setup_name, delta, q_vars,pos_in, pos_out,alpha_target,beta_target,Re_target,Br_target,name_save);


