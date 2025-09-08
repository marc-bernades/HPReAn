function plot_ScaleMotion(y,delta,q_vars,U,S,V,aa,bb,ii,jj,cc,idx,alpha, beta, c_speed, gamma, BF, Dy, label_save, wall_plot)


for n_var = 1:q_vars
    Q{n_var}     = U{aa,bb,ii,jj,cc}(n_var:q_vars:end,idx); % Take all columns for each perturbation row
    F{n_var}     = V{aa,bb,ii,jj,cc}(n_var:q_vars:end,idx); % Take all columns for each perturbation row
end

Sigma(1,cc) = S{aa,bb,ii,jj,cc}(1,1);
Sigma(2,cc) = S{aa,bb,ii,jj,cc}(2,2);
Sigma(3,cc) = S{aa,bb,ii,jj,cc}(3,3);

% Responses normalized with v
Norm  = 1; %max(abs(Q{3}));
rho_r = Q{1}/Norm;
u_r   = Q{2}/Norm;
v_r   = Q{3}/Norm;
w_r   = Q{4}/Norm;
T_r   = Q{5}/Norm;

% Forcing normalized with w
Norm  = 1; %max(abs(F{3}));
rho_f = F{1}/Norm;
u_f   = F{2}/Norm;
v_f   = F{3}/Norm;
w_f   = F{4}/Norm;
T_f   = F{5}/Norm;


%% INPUT-OUTPUT VECTORS TEMPORAL-SPATIAL
% Define x
x      = linspace(0,4*pi*delta,200);  % = Beta X / (2pi)
x_norm = x/delta;   
% Define z
z      = linspace(0,4/3*pi*delta,200); % = Beta Z / (2pi)
z_norm = 0*z/delta;   

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
label_strings       = {'\rho^\prime','u^\prime','v^\prime','w^\prime','T^\prime'};
% label_strings_norm  = {'\rho^\prime / \rho_b','u^\prime / u_b','v^\prime / u_b','w^\prime / u_b','T^\prime / T_b'};

% Output from response
out{1} = (rho_r.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));
out{2} =   (u_r.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));
out{3} =   (v_r.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));
out{4} =   (w_r.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));
out{5} =   (T_r.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));

% Input from forcing
% ph     = -0.0; % IG IN

in{1} = (rho_f.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));
in{2} =   (u_f.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));
in{3} =   (v_f.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));
in{4} =   (w_f.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));
in{5} =   (T_f.*exp(sqrt(-1).*(alpha(aa).*(x_norm + ph_x) + beta(bb).*(z_norm + ph_z) - omega*t)));

% Pseudo-boiling location
y_pb = 1.9 + zeros(size(out{1}(1,:)));


%% RESPONSES OUTPUT
[xx_count,y_count]     = meshgrid(x_norm,y);

% Find y lim at y = 0.5
% [~,y_max] = min(abs(BF.y - 0.5));
% Cond = 1:y_max;
Cond = 1:length(BF.y);

% XY structures for each output
for ii = 1:length(out)
    % XY Contours
    f = figure; hold on; box on; 

    c_min = min(min((real(out{ii}(Cond,:)))));
    c_max = max(max((real(out{ii}(Cond,:)))));

    [c2,h]=contourf(xx_count(Cond,:),y_count(Cond,:),real(out{ii}(Cond,:)),[linspace(c_min,c_max,250)],'HandleVisibility','off');
    c = colorbar('Ticks',[linspace(floor(c_min),floor(c_max)+1,5)]);
    % b_Ticks = [linspace(-3,3,3)]; % LSM
    % b_Ticks = [linspace(-8,8,5)]; % VLSM
    % b_Ticks = [linspace(-4.5,4.5,3)]; % NWSM_tw
%     b_Ticks = [linspace(-3.5,3.5,3)]; % NWSM_bw
    plot(xx_count(1,:), y_pb, 'LineStyle','--','Color','k','LineWidth',1)
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
    xlabel('${x / \delta }$','interpreter','latex','fontsize',14)
    ylabel('${y/\delta}$','interpreter','latex','fontsize',14)
    xticks([0 pi 2*pi 3*pi 4*pi]); xticklabels({'$0$','$\pi$','$2\pi$', '$3\pi$', '$4\pi$'})
    xaxisproperties= get(gca, 'XAxis');
    xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
    if strcmp(wall_plot,'tw')
        ylim([1.8 2]); yticks([1.8 1.9 2.0]);
    else
        ylim([0.0 0.2]); yticks([0.0 0.1 0.2]);
    end
    
    pbaspect([1 0.3 1])
    exportgraphics(gca,strcat('Figures/ScaleMotion_',label_save,'_',label_variable{ii},'_Sigma_',num2str(idx),'_Alpha_',num2str(alpha(aa)),'_Beta_',num2str(beta(bb)),'.jpeg'),'Resolution',300)
    % exportgraphics(gca,strcat('Figures/ResponsePattern_rho_',label_save,'_Sigma_',num2str(idx),'_Alpha_',num2str(alpha(aa)),'_Beta_',num2str(beta(bb)),'.jpeg'),'Resolution',300)

    
end

%% Responses y-normal
f = figure; hold on; box on;

Color_map = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980],...
    [0.4940 0.1840 0.5560], [0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840], ...
    [0.9290 0.6940 0.1250], [1, 0, 0], [0, 0, 1], [0, 1, 0], [1, 1, 0], [0, 0.5, 0.5]};

for ii = 1:length(out)
    var_real = (out{ii}(Cond,1));
    var_name = label_variable{ii};

    subplot(3,1,1)
    plot(y,real(var_real),'Color',Color_map{ii}); hold on;
    set(gca,'fontsize',12)
    set(gca,'linewidth',2.5)
    % xlabel('${y/\delta}$','interpreter','latex','fontsize',14)
    set(gca,'XTick',[])
    ylabel(['$Re(' , 'q^\prime' , ')$'],'interpreter','latex','fontsize',14)


    %pbaspect([1 0.5 1])

    subplot(3,1,2)
    plot(y,imag(var_real),'Color',Color_map{ii}); hold on;
    set(gca,'fontsize',12)
    set(gca,'linewidth',2.5)
    set(gca,'XTick',[])
    % xlabel('${y/\delta}$','interpreter','latex','fontsize',14)
    ylabel(['$Im(' , 'q^\prime' , ')$'],'interpreter','latex','fontsize',14)

    subplot(3,1,3)
    plot(y,abs(var_real),'Color',Color_map{ii}); hold on;
    set(gca,'fontsize',12)
    set(gca,'linewidth',2.5)
    xlabel('${y/\delta}$','interpreter','latex','fontsize',14)
    ylabel(['$\vert ' , 'q^\prime' , ' \vert$'],'interpreter','latex','fontsize',14)

    c_legend{ii} = ['$',label_strings{ii},'$'];
end


legend(c_legend,'interpreter','latex','NumColumns',3,'fontsize',12,'Location','northwest',Box='off')

exportgraphics(gca,strcat('Figures/ScaleMotion_Y_',label_save,'_Sigma_',num2str(idx),'_Alpha_',num2str(alpha(aa)),'_Beta_',num2str(beta(bb)),'.jpeg'),'Resolution',300)



% %% Superimpoisition
% u_base_out = repmat(BF.u_star,1,length(xx_count(2,:)));
% 
% % f = figure
% % c_min = min(min((real(u_base_out))));
% % c_max = max(max((real(u_base_out))));
% % [c1,h1]=contourf(xx_count,y_count,real(u_base_out),[linspace(c_min,c_max,5)],'HandleVisibility','off');
% 
% 
% y_target    = 1.90;
% y_target    = 0.22;
% 
% [~,idx_y]   = (min(abs(y_count(:,2) - (y_target)))); % y_count(idx_y,2); % Position pseudo
% xi          = 0.1;
% idx_y_off   = 2; %2 LSM tw, -5 VSM / bw; % offset superimposition range
% u_super_out = u_base_out;
% u_super_out(idx_y-idx_y_off:end,:) = u_base_out(idx_y-idx_y_off:end,:) + xi*real(u_out(idx_y-idx_y_off:end,:));
% % u_super_out(1:idx_y-idx_y_off,:) = u_base_out(1:idx_y-idx_y_off,:) + xi*real(u_out(1:idx_y-idx_y_off,:));
% 
% f = figure
% c_min = min(min((real(u_super_out))));
% c_max = max(max((real(u_super_out))));
% [c2,h2]=contourf(xx_count,y_count,real(u_super_out),[linspace(c_min,c_max,8)],'HandleVisibility','off');
% c2 = colorbar('Ticks',[linspace(c_min,c_max,5)]);
% % set(h2, 'edgecolor','none');
% colorbar
% colormap('jet')
% cbh = findall(f, 'Type', 'ColorBar');
% cTH = get(cbh,'Title');
% set(cTH,'String',['$','u_0 + u^\prime','$'],'Interpreter','latex','fontsize',14);
% pbaspect([1 1.0 1])
% set(gca,'fontsize',12)
% set(gca,'linewidth',1.5)
% xlabel('${x / \delta }$','interpreter','latex','fontsize',14)
% ylabel('${y/\delta}$','interpreter','latex','fontsize',14)
% xticks([0 pi 2*pi 3*pi 4*pi]); xticklabels({'$0$','$\pi$','$2\pi$', '$3\pi$', '$4\pi$'})
% xaxisproperties= get(gca, 'XAxis');
% xaxisproperties.TickLabelInterpreter = 'latex'; % latex for x-axis
% pbaspect([1 0.3 1])
% 
% % exportgraphics(gca,strcat('Figures/ScaleMotion_Superimposed_',label_save,'_Xi_',num2str(xi),'_Sigma_',num2str(idx),'_Alpha_',num2str(alpha(aa)),'_Beta_',num2str(beta(bb)),'.jpeg'),'Resolution',300)




% f = figure;
% u_base_y = BF.u(idx_y)./BF.norm.u;
% plot(xx_count(idx_y,:),u_base_y + xi*u_out(idx_y,:))
% set(gca,'fontsize',12)
% set(gca,'linewidth',1.5)
% xlabel('${x/\delta}$','interpreter','latex','fontsize',14)
% ylabel('${u / u_b }$','interpreter','latex','fontsize',14)
% pbaspect([1 0.3 1])








