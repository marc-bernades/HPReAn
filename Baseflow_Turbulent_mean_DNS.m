function [BF, Dy] = Baseflow_Turbulent_mean_DNS(y,D,delta,bScaling)

%% Load DNS dataset
% DNS 128^3
source_path  = '/home/marc/Documents/Doctorat/Solver_and_postprocess/RHEA/Paper7/';
filename     = '3d_high_pressure_turbulent_channel_flow_28500000.h5';

DNS{1}       = Load_dataset(source_path, filename);
DNS{1}.save  = 'DNS';
DNS{1}.label = '$ DNS $';
delta_h      = 100*1E-6;


%% Loop DNS files
for ii = 1:length(DNS)
    % y-discrete
    BF{ii}.y     = DNS{ii}.y(2,:,2);
    % Spatial average fields
    BF{ii}.u     = Spatial_avg_XZ_var(DNS{ii}.avg_u,DNS{ii}.num_points_x,DNS{ii}.num_points_y,DNS{ii}.num_points_z, 0);
    BF{ii}.rho   = Spatial_avg_XZ_var(DNS{ii}.avg_rho,DNS{ii}.num_points_x,DNS{ii}.num_points_y,DNS{ii}.num_points_z, 0);
    BF{ii}.T     = Spatial_avg_XZ_var(DNS{ii}.avg_T,DNS{ii}.num_points_x,DNS{ii}.num_points_y,DNS{ii}.num_points_z, 0);
    BF{ii}.c_p   = Spatial_avg_XZ_var(DNS{ii}.avg_c_p,DNS{ii}.num_points_x,DNS{ii}.num_points_y,DNS{ii}.num_points_z, 0);
    BF{ii}.c_v   = Spatial_avg_XZ_var(DNS{ii}.avg_c_v,DNS{ii}.num_points_x,DNS{ii}.num_points_y,DNS{ii}.num_points_z, 0);
    BF{ii}.mu    = Spatial_avg_XZ_var(DNS{ii}.avg_mu,DNS{ii}.num_points_x,DNS{ii}.num_points_y,DNS{ii}.num_points_z, 0);
    BF{ii}.kappa = Spatial_avg_XZ_var(DNS{ii}.avg_kappa,DNS{ii}.num_points_x,DNS{ii}.num_points_y,DNS{ii}.num_points_z, 0);
    BF{ii}.P_0   = Spatial_avg_XZ_var(DNS{ii}.avg_P,DNS{ii}.num_points_x,DNS{ii}.num_points_y,DNS{ii}.num_points_z, 0);
    BF{ii}.E     = Spatial_avg_XZ_var(DNS{ii}.avg_E,DNS{ii}.num_points_x,DNS{ii}.num_points_y,DNS{ii}.num_points_z, 0);
    DNS{ii}.avg_ke   = 1/2.*(DNS{ii}.u.^2 + DNS{ii}.v.^2 + DNS{ii}.w.^2);
    DNS{ii}.avg_e    = DNS{ii}.avg_E - DNS{ii}.avg_ke;
    BF{ii}.e     = Spatial_avg_XZ_var(DNS{ii}.avg_e,DNS{ii}.num_points_x,DNS{ii}.num_points_y,DNS{ii}.num_points_z, 0);
    BF{ii}.sos   = Spatial_avg_XZ_var(DNS{ii}.avg_sos,DNS{ii}.num_points_x,DNS{ii}.num_points_y,DNS{ii}.num_points_z, 0);
    BF{ii}.ke    = Spatial_avg_XZ_var(DNS{ii}.avg_ke,DNS{ii}.num_points_x,DNS{ii}.num_points_y,DNS{ii}.num_points_z, 0);

    BF{ii}.lambda = -2/3*BF{ii}.mu;

    % Calculate bulk values, cell size and deltas
    [BF{ii}.rho_b, BF{ii}.mu_b, BF{ii}.P_b, BF{ii}.T_b, BF{ii}.u_b, BF{ii}.kappa_b, BF{ii}.c_p_b, ...
        BF{ii}.c_v_b, BF{ii}.E_b, BF{ii}.e_b, BF{ii}.sos_b, ...
        BF{ii}.Re_b, BF{ii}.Pr_b, BF{ii}.Ec_b, BF{ii}.Ma_b] = Calculate_bulk_baseflow(DNS{ii}.avg_rho, DNS{ii}.avg_mu, DNS{ii}.avg_P, DNS{ii}.avg_T, DNS{ii}.avg_u, DNS{ii}.avg_kappa, DNS{ii}.avg_c_p, DNS{ii}.avg_c_v, DNS{ii}.avg_E, DNS{ii}.avg_e, DNS{ii}.avg_sos, DNS{ii}.num_points_x, DNS{ii}.num_points_y, DNS{ii}.num_points_z, delta_h, DNS{ii}.dx, DNS{ii}.dy, DNS{ii}.dz);

    % Calculate wall values (based on time-averaged last dataset)
    [BF{ii}.rho_bw, BF{ii}.mu_bw, BF{ii}.T_tau_bw, BF{ii}.T_bw, BF{ii}.u_tau_bw,BF{ii}.rho_tw, BF{ii}.mu_tw, BF{ii}.T_tau_tw, BF{ii}.T_tw, BF{ii}.u_tau_tw, ...
        BF{ii}.c_p_bw, BF{ii}.kappa_bw, BF{ii}.c_p_tw, BF{ii}.kappa_tw,...
        BF{ii}.Re_tau_bw, BF{ii}.Cf_bw, BF{ii}.Pr_bw, BF{ii}.Nu_bw, BF{ii}.St_bw, BF{ii}.Re_tau_tw, BF{ii}.Cf_tw, BF{ii}.Pr_tw, BF{ii}.Nu_tw, BF{ii}.St_tw] = Calculate_wall_metrics(DNS{ii}.x,DNS{ii}.y,DNS{ii}.z,DNS{ii}.avg_rho,DNS{ii}.avg_mu,DNS{ii}.avg_T,DNS{ii}.avg_c_p,DNS{ii}.avg_kappa,DNS{ii}.avg_u,...
        BF{ii}.rho_b,BF{ii}.T_b,BF{ii}.u_b,DNS{ii}.num_points_x,DNS{ii}.num_points_y,DNS{ii}.num_points_z,delta_h);
    
    % Normalize wall-units (based on time-averaged last dataset)
    [BF{ii}.y_plus_bw,BF{ii}.u_plus_bw,BF{ii}.T_plus_bw,BF{ii}.y_plus_tw,BF{ii}.u_plus_tw,BF{ii}.T_plus_tw] = Transform_WallUnits(DNS{ii}.y,DNS{ii}.avg_u,DNS{ii}.avg_T,BF{ii}.u_tau_bw,BF{ii}.rho_bw,BF{ii}.mu_bw,BF{ii}.T_bw,BF{ii}.T_tau_bw,BF{ii}.u_tau_tw,BF{ii}.rho_tw,BF{ii}.mu_tw,BF{ii}.T_tw,BF{ii}.T_tau_tw,DNS{ii}.num_points_x,DNS{ii}.num_points_y,DNS{ii}.num_points_z,delta_h);


end

% Interpolate y-direction to Chebyshev
struct_field_names = fieldnames(BF{1});
y_DNS = BF{1}.y/delta_h;
for kk = 1:length(struct_field_names)
    try
        BF{1}.(struct_field_names{kk}) = interpn(y_DNS,BF{1}.(struct_field_names{kk}),y')';
    catch
        aa = 1;
    end

end

BF{1}.y       = y; % To accomodate units
BF{1}.F_hat_0 = 0;
BF{1}.delta_h = delta_h;

% Derivative operator
Dy{1}.D          = D;
Dy{1}.du_dy      = D*BF{1}.u;      %CentralDerivative_d1_2ndOrder_1D(u,y);
Dy{1}.dT_dy      = D*BF{1}.T;      %CentralDerivative_d1_2ndOrder_1D(T,y);
Dy{1}.dmu_dy     = D*BF{1}.mu;     %CentralDerivative_d1_2ndOrder_1D(mu,y);
Dy{1}.dkappa_dy  = D*BF{1}.kappa;  %CentralDerivative_d1_2ndOrder_1D(kappa,y);
Dy{1}.drho_dy    = D*BF{1}.rho;    %CentralDerivative_d1_2ndOrder_1D(rho,y);
Dy{1}.dlambda_dy = D*BF{1}.lambda; %CentralDerivative_d1_2ndOrder_1D(lambda,y);
Dy{1}.de_dy      = D*BF{1}.e;      %CentralDerivative_d1_2ndOrder_1D(e,y);
Dy{1}.du2_dy2    = D*D*BF{1}.u;    %CentralDerivative_d2_2ndOrder_1D(u,y);
Dy{1}.dT2_dy2    = D*D*BF{1}.T;    %CentralDerivative_d2_2ndOrder_1D(T,y);

% Normalization
BF{1}.norm   = Reference_Scalings(BF{1}.u_b,BF{1}.T_b,BF{1}.rho_b,BF{1}.mu_b,BF{1}.kappa_b,BF{1}.c_p_b,BF{1}.c_v_b,BF{1}.E_b,BF{1}.e_b,BF{1}.P_b,...
    BF{1}.u,BF{1}.T,BF{1}.rho,BF{1}.mu,BF{1}.kappa,BF{1}.c_p,BF{1}.c_v,BF{1}.lambda,Dy,delta,bScaling, BF{1}.u_b);

BF{1}.u_star     = BF{1}.u./BF{1}.norm.u;
BF{1}.T_star     = BF{1}.T./BF{1}.norm.T;
BF{1}.kappa_star = BF{1}.kappa./BF{1}.norm.kappa;
BF{1}.mu_star    = BF{1}.mu./BF{1}.norm.mu;

% Add non_b quantities
BF{1}.Re = BF{1}.Re_b;
BF{1}.Pr = BF{1}.Pr_b;
BF{1}.Ec = BF{1}.Ec_b;
BF{1}.Ma = BF{1}.Ma_b;

end