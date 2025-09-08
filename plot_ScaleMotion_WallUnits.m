function plot_ScaleMotion_WallUnits(y,delta,q_vars,U,S,V,aa,bb,ii,jj,cc,idx,alpha, beta, c_speed, gamma, BF, Dy, label_save, wall_plot)


for n_var = 1:q_vars
    Q{n_var}     = U{aa,bb,ii,jj,cc}(n_var:q_vars:end,idx); % Take all columns for each perturbation row
    F{n_var}     = V{aa,bb,ii,jj,cc}(n_var:q_vars:end,idx); % Take all columns for each perturbation row
end

Sigma(1,cc) = S{aa,bb,ii,jj,cc}(1,1);
Sigma(2,cc) = S{aa,bb,ii,jj,cc}(2,2);
Sigma(3,cc) = S{aa,bb,ii,jj,cc}(3,3);

% Convert to u+
Norm{1}  = BF.(strcat('rho_', wall_plot))/BF.norm.rho;
Norm{2}  = BF.(strcat('u_tau_', wall_plot))/BF.norm.u;
% Q{5}     = Q{5}*BF.norm.T;
% T_r      = abs((Q{5} - BF.(strcat('T_', wall_plot)))/BF.(strcat('T_tau_', wall_plot)));
Norm{3}  = BF.(strcat('T_tau_', wall_plot))/BF.norm.T;

% Responses normalized with v
rho_r = Q{1}/Norm{1};
u_r   = Q{2}/Norm{2};
v_r   = Q{3}/Norm{2};
w_r   = Q{4}/Norm{2};
T_r   = Q{5}/Norm{3};




%% INPUT-OUTPUT VECTORS TEMPORAL-SPATIAL
% Define x
x      = linspace(0,4*pi*delta,200);  % = Beta X / (2pi)
x_norm = x;   
% Define z
z      = linspace(0,4/3*pi*delta,200); % = Beta Z / (2pi)
z_norm = 0*z;   

% [x_norm, z_norm] = meshgrid(x_norm,z_norm);

% Phase shift
ph_x     = 0.00; %HP
ph_z     = 0.00; %HP

% time
omega = c_speed(cc)*alpha(aa);
t = -gamma/(omega);

% Adjust wiggles
% u_r = smooth(u_r(1:1:end),3);

% Convert to x,y,z,t space
label_variable = {'rho','u','v','w','T'};
label_strings  = {'\rho^\prime/\rho_w','u^{\prime +}','v^{\prime +}','w^{\prime +}','T^{\prime +}'};

% Output from response
out{1} = (rho_r.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));
out{2} =   (u_r.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));
out{3} =   (v_r.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));
out{4} =   (w_r.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));
out{5} =   (T_r.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));


% Pseudo-boiling location
y_pb = 1.9 + zeros(size(out{1}(1,:)));
y_pb_plus = BF.delta_h*(2*delta - y_pb)*(BF.u_tau_tw/(BF.mu_tw/BF.rho_tw));

%% RESPONSES OUTPUT
y_2_yplus = (BF.(strcat('u_tau_', wall_plot))/(BF.(strcat('mu_', wall_plot))/BF.(strcat('rho_', wall_plot))))*BF.delta_h;

if strcmp(wall_plot,'tw')
    y         = (max(y) - y)*y_2_yplus;
else
    y = y*y_2_yplus;
end
x_norm    = x_norm*y_2_yplus;

[xx_count,y_count]     = meshgrid(x_norm,y);

% Find y lim at y+ = 100
% [~,y_plus_100] = min(abs(y - 100));
Cond = 1:length(y);

% XY structures for each output
for ii = 1:length(out)
    % XY Contours
    f = figure; hold on; box on;

    c_min = min(min((real(out{ii}(Cond,:)))));
    c_max = max(max((real(out{ii}(Cond,:)))));

    [c2,h]=contourf(xx_count(Cond,:),y_count(Cond,:),real(out{ii}(Cond,:)),[linspace(c_min,c_max,250)],'HandleVisibility','off');
    c = colorbar('Ticks',[linspace(floor(c_min),floor(c_max)+1,5)]);
    set(gca,'YScale', 'log')

    if strcmp(wall_plot,'tw')
        plot(xx_count(1,:), y_pb_plus, 'LineStyle','--','Color','k','LineWidth',1)
    end
    % c.Location = 'northoutside';
    %     clim([b_Ticks(1) b_Ticks(end)]);
    %     c = colorbar('Ticks',b_Ticks);
    set(h, 'edgecolor','none');
    % colorbar
    colormap(redblue(250)) % colormap('jet')
    cbh = findall(f, 'Type', 'ColorBar');
    cTH = get(cbh,'Title');
    set(cTH,'String',['$',label_strings{ii},'$'],'Interpreter','latex','fontsize',14);
    pbaspect([1 1.0 1])
    set(gca,'fontsize',12)
    set(gca,'linewidth',1.5)
    xlabel('${x^+}$','interpreter','latex','fontsize',14)
    ylabel('${y^+}$','interpreter','latex','fontsize',14)
    xlim([0 1000])
    xticks([0 250 500 750 1000]); xticklabels({'$0$','$250$','$500$','$750$','$1000$'})
    xaxisproperties= get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
    ylim([0.1 100])
    yticks([0.1 1 10 100]); yticklabels({'$10^{-1}$', '$10^0$','$10^1$','$10^2$'})
    yaxisproperties= get(gca, 'YAxis');
    yaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
    pbaspect([1 0.3 1])

    exportgraphics(gca,strcat('Figures/ScaleMotion_WallUnits_',label_save,'_',label_variable{ii},'_Sigma_',num2str(idx),'_Alpha_',num2str(alpha(aa)),'_Beta_',num2str(beta(bb)),'.jpeg'),'Resolution',300)


end


% Streamwise velocity
% f = figure;
% 
% c_min = min(min((real(u_out(Cond,:)))));
% c_max = max(max((real(u_out(Cond,:)))));
% 
% [c2,h]=contourf(xx_count(Cond,:),y_count(Cond,:),real(u_out(Cond,:)),[linspace(c_min,c_max,250)],'HandleVisibility','off');
% % set(gca,'XScale', 'log')
% set(gca,'YScale', 'log')
% 
% b_Ticks = [linspace(-50,50,3)]; % LSM 
% % c.Location = 'northoutside';
% clim([b_Ticks(1) b_Ticks(end)]);
% c = colorbar('Ticks',b_Ticks);
% set(h, 'edgecolor','none');
% % colorbar
% colormap(redblue(250)) % colormap('jet')
% cbh = findall(f, 'Type', 'ColorBar');
% cTH = get(cbh,'Title');
% set(cTH,'String',['$','{u^{\prime \thinspace +}}','$'],'Interpreter','latex','fontsize',14);
% pbaspect([1 1.0 1])
% set(gca,'fontsize',12)
% set(gca,'linewidth',1.5)
% xlabel('${x^+}$','interpreter','latex','fontsize',14)
% ylabel('${y^+}$','interpreter','latex','fontsize',14)
% xlim([0 1000])
% xticks([0 250 500 750 1000]); xticklabels({'$0$','$250$','$500$','$750$','$1000$'})
% xaxisproperties= get(gca, 'XAxis');
% xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
% ylim([0.1 100])
% yticks([0.1 1 10 100]); yticklabels({'$10^{-1}$', '$10^0$','$10^1$','$10^2$'})
% yaxisproperties= get(gca, 'YAxis');
% yaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
% pbaspect([1 0.3 1])
% exportgraphics(gca,strcat('Figures/ScaleMotion_WallUnits_',label_save,'_Sigma_',num2str(idx),'_Alpha_',num2str(alpha(aa)),'_Beta_',num2str(beta(bb)),'.jpeg'),'Resolution',300)
% % exportgraphics(gca,strcat('Figures/ResponsePattern_rho_',label_save,'_Sigma_',num2str(idx),'_Alpha_',num2str(alpha(aa)),'_Beta_',num2str(beta(bb)),'.jpeg'),'Resolution',300)
% 
% 
% 
% f = figure;
% plot(y(1:y_plus_100),u_r(1:y_plus_100))
% set(gca,'fontsize',12)
% set(gca,'linewidth',1.5)
% xlabel('${y/\delta}$','interpreter','latex','fontsize',14)
% ylabel('${u / u_b }$','interpreter','latex','fontsize',14)
% pbaspect([1 0.5 1])
% exportgraphics(gca,strcat('Figures/ScaleMotion',label_save,'u_r_Sigma_',num2str(idx),'_Alpha_',num2str(alpha(aa)),'_Beta_',num2str(beta(bb)),'.jpeg'),'Resolution',300)






