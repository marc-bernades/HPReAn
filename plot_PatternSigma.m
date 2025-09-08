function plot_PatternSigma(y,delta,q_vars,U,S,V,aa,bb,ii,jj,cc,idx,alpha, beta, c_speed, label_save)


for n_var = 1:q_vars
    Q{n_var}     = U{aa,bb,ii,jj,cc}(n_var:q_vars:end,idx); % Take all columns for each perturbation row
    F{n_var}     = V{aa,bb,ii,jj,cc}(n_var:q_vars:end,idx); % Take all columns for each perturbation row
end

Sigma(1,cc) = S{aa,bb,ii,jj,cc}(1,1);
Sigma(2,cc) = S{aa,bb,ii,jj,cc}(2,2);
Sigma(3,cc) = S{aa,bb,ii,jj,cc}(3,3);

% Responses normalized with v
Norm  = max(abs(Q{3}));
rho_r = Q{1}/Norm;
u_r   = Q{2}/Norm;
v_r   = Q{3}/Norm;
w_r   = Q{4}/Norm;
T_r   = Q{5}/Norm;

% Forcing normalized with w
Norm  = max(abs(F{3}));
rho_f = F{1}/Norm;
u_f   = F{2}/Norm;
v_f   = F{3}/Norm;
w_f   = F{4}/Norm;
T_f   = F{5}/Norm;


%% INPUT-OUTPUT VECTORS TEMPORAL-SPATIAL
% Define z
z      = linspace(0,1,200); % = Beta Z / (2pi)
z_norm = z*2*pi./beta(bb);   

% Phase shift
ph     = 0.00; %HP
% ph     = 0.50; %IG
% ph     = 0.2; %Turbulent


% time
t = 0;

% Convert to x,y,z,t space
% Output from response
rho_out = (rho_r.*exp(sqrt(-1).*beta(bb).*(z_norm + ph) + sqrt(-1)*c_speed(cc)*alpha(aa)*t));
u_out   =   (u_r.*exp(sqrt(-1).*beta(bb).*(z_norm + ph) + sqrt(-1)*c_speed(cc)*alpha(aa)*t));
v_out   =   (v_r.*exp(sqrt(-1).*beta(bb).*(z_norm + ph) + sqrt(-1)*c_speed(cc)*alpha(aa)*t));
w_out   =   (w_r.*exp(sqrt(-1).*beta(bb).*(z_norm + ph) + sqrt(-1)*c_speed(cc)*alpha(aa)*t));
T_out   =   (T_r.*exp(sqrt(-1).*beta(bb).*(z_norm + ph) + sqrt(-1)*c_speed(cc)*alpha(aa)*t));
vw_out  = sqrt((v_out).^2 + real(w_out).^2);

% Input from forcing
% ph     = -0.0; % IG IN

rho_in = (rho_f.*exp(sqrt(-1).*beta(bb).*(z_norm + ph) + sqrt(-1)*c_speed(cc)*alpha(aa)*t));
u_in   =   (u_f.*exp(sqrt(-1).*beta(bb).*(z_norm + ph) + sqrt(-1)*c_speed(cc)*alpha(aa)*t));
v_in   =   (v_f.*exp(sqrt(-1).*beta(bb).*(z_norm + ph) + sqrt(-1)*c_speed(cc)*alpha(aa)*t));
w_in   =   (w_f.*exp(sqrt(-1).*beta(bb).*(z_norm + ph) + sqrt(-1)*c_speed(cc)*alpha(aa)*t));
T_in   =   (T_f.*exp(sqrt(-1).*beta(bb).*(z_norm + ph) + sqrt(-1)*c_speed(cc)*alpha(aa)*t));
vw_in  = sqrt((v_in).^2 + real(w_in).^2);


%% FORCING INPUT
i_plot     = 8;

figure
ax = axes('Color', [1 1 1]); 
hold(ax,'on')
box(ax,'on')
% Define the colormap for the quiver arrows.
% cmap can have any number of rows.
cmap        = jet(255); 
ax.Colormap = cmap; 
% Assign colors based on magnitude of vectors
vectorMagnitude = abs(vw_in); 
% Scale magnitudes to rows of colormap
vecMagNorm = (vectorMagnitude-min(vectorMagnitude))./range(vectorMagnitude);
vecColorIdx = round(vecMagNorm * (size(cmap,1)-1)) + 1; 
% Plot the quiver data
for idx_z = 1:i_plot:numel(z_norm)
    for idx_y = 1:i_plot:numel(y)
    quiver(ax, z_norm(idx_z),y(idx_y),real(w_in(idx_y,idx_z)),real(v_in(idx_y,idx_z)), ...
        .15,'Color', cmap(vecColorIdx(idx_y,idx_z),:), 'LineWidth',1.0)
    end
end

% Set properties for the main axes
% axis equal
xlim(ax, [0 max(z_norm)])%xlim(ax, [0 beta(bb)])
ylim(ax, [0 2])
xlabel('${\lambda_z \thinspace / \thinspace (\delta \thinspace k_z)}$','interpreter','latex','fontsize',14)
ylabel('${y/\delta}$','interpreter','latex','fontsize',14)
xticks([0 pi/2 pi]); xticklabels({'$0$','$\pi/2$','$\pi$'})
xaxisproperties= get(ax, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis


% Add colorbar
cb = colorbar(ax); 
% Set colorbar range
caxis(ax, [floor(vectorMagnitude(1)), ceil(vectorMagnitude(2))])
cb.Ticks = linspace(floor(vectorMagnitude(1)),ceil(vectorMagnitude(2)),5);
% Label the colorbars
% ylabel(cb,'Vector magnitude')

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)


cTH = get(cb,'Title');
set(cTH,'String',['$','\vert \textbf{u}^\prime \vert','$'],'Interpreter','latex','fontsize',14);
pbaspect([1 1.0 1])

exportgraphics(gca,strcat('Figures/Forcing_vw_Pattern_',label_save,'_Sigma_',num2str(idx),'_Alpha_',num2str(alpha(aa)),'_Beta_',num2str(beta(bb)),'.jpeg'),'Resolution',300)

%% RESPONSES OUTPUT
[zz_count,y_count]     = meshgrid(z_norm,y);
% Density
f = figure;
c_min = min(min((real(rho_out))));
c_max = max(max((real(rho_out))));

[c,h]=contourf(zz_count,y_count,real(rho_out),[linspace(c_min,c_max,250)],'HandleVisibility','off');
c = colorbar('Ticks',[linspace(c_min,c_max,5)]);
%     clim([b_Ticks(1) b_Ticks(end)]);
%     c = colorbar('Ticks',b_Ticks);
set(h, 'edgecolor','none');
colorbar
colormap(redblue(250)) % colormap('jet')
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','\rho^\prime','$'],'Interpreter','latex','fontsize',14);
pbaspect([1 1.0 1])
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
xlabel('${\lambda_z \thinspace / \thinspace (\delta \thinspace k_z)}$','interpreter','latex','fontsize',14)
ylabel('${y/\delta}$','interpreter','latex','fontsize',14)
xticks([0 pi/2 pi]); xticklabels({'$0$','$\pi/2$','$\pi$'})
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis


exportgraphics(gca,strcat('Figures/ResponsePattern_rho_',label_save,'_Sigma_',num2str(idx),'_Alpha_',num2str(alpha(aa)),'_Beta_',num2str(beta(bb)),'.jpeg'),'Resolution',300)


% Streamwise
f = figure;
c_min = min(min((real(u_out))));
c_max = max(max((real(u_out))));

[c,h]=contourf(zz_count,y_count,real(u_out),[linspace(c_min,c_max,250)],'HandleVisibility','off');
c = colorbar('Ticks',[linspace(c_min,c_max,5)]);
%     clim([b_Ticks(1) b_Ticks(end)]);
%     c = colorbar('Ticks',b_Ticks);
set(h, 'edgecolor','none');
colorbar
colormap(redblue(250)) % colormap('jet')
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','u^\prime','$'],'Interpreter','latex','fontsize',14);
pbaspect([1 1.0 1])
set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
xlabel('${\lambda_z \thinspace / \thinspace (\delta \thinspace k_z)}$','interpreter','latex','fontsize',14)
ylabel('${y/\delta}$','interpreter','latex','fontsize',14)
xticks([0 pi/2 pi]); xticklabels({'$0$','$\pi/2$','$\pi$'})
xaxisproperties= get(gca, 'XAxis');
xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis


exportgraphics(gca,strcat('Figures/ResponsePattern_u_',label_save,'_Sigma_',num2str(idx),'_Alpha_',num2str(alpha(aa)),'_Beta_',num2str(beta(bb)),'.jpeg'),'Resolution',300)


