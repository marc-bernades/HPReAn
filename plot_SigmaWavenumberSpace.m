function plot_SigmaWavenumberSpace(alpha, beta, N_target, n_LST_Sweep, omega, Sigma_1,  ...
    Sigma_1_min, Sigma_1_max, n_steps,name_file_load, b_gain_log, b_axis_labels_log, b_Ticks, varargin)

% Contourplot interpolation
f = figure; hold on; grid on; box on;

alpha_q          = linspace(min(alpha),max(alpha),500);
beta_q           = linspace(min(beta),max(beta),500);
if length(b_gain_log) > 1
    alpha_q = logspace(log10(min(alpha)),log10(max(alpha)),500);
    beta_q  = logspace(log10(min(beta)),log10(max(beta)),500);
end
[beta,alpha]     = meshgrid(beta,alpha);
[beta_q,alpha_q] = meshgrid(beta_q,alpha_q);
Sigma_1_q          = griddata(alpha,beta,Sigma_1,alpha_q,beta_q,'linear'); %G_max;
Sigma_1_q(Sigma_1_q<Sigma_1_min) = Sigma_1_min;


% Contourplot
% [c,h]=contourf(alpha_q,beta_q,Sigma_1_q,[linspace(min(min(Sigma_1_q)),round(max(max(Sigma_1_q)),1),500)],'HandleVisibility','off');

if length(b_gain_log) == 1
    [c,h]=contourf(alpha_q,beta_q,Sigma_1_q,[linspace(Sigma_1_min,Sigma_1_max,20)],'HandleVisibility','off');
    c = colorbar('Ticks',[linspace(Sigma_1_min,Sigma_1_max,n_steps)]);
    clim([b_Ticks(1) b_Ticks(end)]);
    c = colorbar('Ticks',b_Ticks);


else % Logs
    log_min = log10(Sigma_1_min); if isinf(log_min) log_min = -2; end
    log_max = round(log10(Sigma_1_max));
    [c,h]=contourf(alpha,beta,(Sigma_1),[logspace(log_min,log_max,20)],'HandleVisibility','off');
    set(gca,'ColorScale','log')
    c = colorbar('Ticks',[logspace(b_gain_log(1),b_gain_log(2),diff(b_gain_log)+1)]);
    clim([10^b_gain_log(1) 10^b_gain_log(2)])

end
set(h, 'edgecolor','none');
if length(b_axis_labels_log) > 1
    set(gca,'XScale', 'log')
    set(gca,'YScale', 'log')
end

% colorbar;
colormap('jet')
% c.Ruler.Exponent        = floor(log10(GR_max));                   % set the desired exponent
% c.Ruler.TickLabelFormat = '%0.1f';       % fix up ugly default %g formatting
cbh = findall(f, 'Type', 'ColorBar');
cTH = get(cbh,'Title');
set(cTH,'String',['$','\sigma_{1}','$'],'Interpreter','latex','fontsize',14);
% set(cTH,'String',['$','log_{10}(\sigma_{1})','$'],'Interpreter','latex','fontsize',14);
% set(cTH,'String',['$','log_{10}(\sigma_{1})','$'],'Interpreter','latex','fontsize',14);
% set(cTH,'String',['$','g','$'],'Interpreter','latex','fontsize',14);

pbaspect([1 1.5 1])

set(gca,'fontsize',12)
set(gca,'linewidth',1.5)
xlabel('${k_x}$','interpreter','latex')
ylabel('${k_z}$','interpreter','latex')
% ylim([0 5]); yticks([0.0 1.0 2.0 3.0 4.0 5.0]);
set(cTH,'Units','normalized','Position',[0.5 1.010 0])

% if size(alpha,1)>1
%     xlim([min(alpha(:,1)) max(alpha(:,1))]); %xticklabels({'2','4','6','8','10'})
% else
%     xlim([min(alpha) max(alpha)]); %xticklabels({'2','4','6','8','10'})
% end

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)

exportgraphics(gca,strcat('Figures/Sigma_Map_',name_file_load,'_c_',num2str(omega),'_Br_',num2str(N_target),'_Re_',num2str(n_LST_Sweep),'.jpeg'),'Resolution',300)



% semilogx(alpha(:,1),Sigma_1(:,1))
% semilogx(beta(1,:),Sigma_1(1,:))


end