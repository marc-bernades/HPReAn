function [rho_b, mu_b, P_b, T_b, u_b, kappa_b, c_p_b, c_v_b, E_b, e_b, sos_b, ...
    Re_b, Pr_b, Ec_b, Ma_b] = Calculate_bulk_baseflow(avg_rho, avg_mu, avg_P, avg_T, avg_u, avg_kappa, avg_c_p, avg_c_v, avg_E, avg_e, avg_sos, num_points_x, num_points_y, num_points_z, delta_h, dx, dy, dz)

delta = zeros(size(avg_rho)); 
volume_total = 0;
rho_b = 0;
mu_b  = 0;
P_b   = 0;
T_b   = 0;
u_b   = 0;
kappa_b = 0;
c_p_b = 0;
c_v_b = 0;
E_b   = 0;
e_b   = 0;
sos_b = 0;


for ii = 2:num_points_x-1
    for jj = 2:num_points_y-1
        for kk = 2:num_points_z-1


            volume          = dx(ii,jj,kk)*dy(ii,jj,kk)*dz(ii,jj,kk);

            delta(ii,jj,kk) = volume.^(1/3);

            % Bulk values
            volume_total    = volume_total + volume;
            rho_b = rho_b + avg_rho(ii,jj,kk)*volume;
            mu_b  = mu_b  + avg_mu(ii,jj,kk)*volume;
            P_b   = P_b   + avg_P(ii,jj,kk)*volume;
            T_b   = T_b   + avg_T(ii,jj,kk)*volume;
            u_b   = u_b   + avg_u(ii,jj,kk)*volume;
            kappa_b   = kappa_b   + avg_kappa(ii,jj,kk)*volume;
            c_p_b     = c_p_b   + avg_c_p(ii,jj,kk)*volume;
            c_v_b     = c_v_b   + avg_c_v(ii,jj,kk)*volume;
            E_b       = E_b   + avg_E(ii,jj,kk)*volume;
            e_b       = e_b   + avg_e(ii,jj,kk)*volume;
            sos_b     = sos_b   + avg_sos(ii,jj,kk)*volume;


        end
    end
end


% Normalize bulk values
rho_b = rho_b/volume_total;
mu_b  = mu_b/volume_total;
P_b   = P_b/volume_total;
T_b   = T_b/volume_total;
u_b   = u_b/volume_total;
kappa_b   = kappa_b/volume_total;
c_p_b = c_p_b/volume_total;
c_v_b = c_v_b/volume_total;
E_b   = E_b/volume_total;
e_b   = e_b/volume_total;
sos_b = sos_b/volume_total;

% Dimensionless numbers
Re_b = u_b*(delta_h)*rho_b/mu_b;
Pr_b = mu_b.*c_p_b/kappa_b;
Ec_b = u_b^2./(c_p_b*T_b);
Ma_b = u_b/sos_b;



end