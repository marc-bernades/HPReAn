%% Baseflow comparison
% HP Poiseuille
% IG Poiseuille
% DNS

% Poiseuille
File_sweep{1} =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_HP';
File_sweep{2} =  'N2_Tbw_94.6425_Ttw_189.285_Pb_6791600_PrEc_bulk_IG';

for ff = 1:length(File_sweep)
    % Unpack structure and load variables
    Results{ff} = load(strcat('Results/', File_sweep{ff}, '.mat'));
end

% DNS
Main_DNS;

Color_map = {[0 0.4470 0.7410],[0.4660 0.6740 0.1880],[0.8500 0.3250 0.0980],...
    [0.4940 0.1840 0.5560], [0.3010 0.7450 0.9330],[0.6350 0.0780 0.1840], ...
    [0.9290 0.6940 0.1250], [1, 0, 0], [0, 0, 1], [0, 1, 0], [1, 1, 0], [0, 0.5, 0.5]};

% Streamwise plot
f = figure
hold on; grid on; box on;
plot(Results{1}.Data.y,Results{1}.Data.BF{1}.u/Results{1}.Data.BF{1}.norm.u,'LineWidth',2,'LineStyle','-','color',Color_map{1});
plot(Results{2}.Data.y,Results{2}.Data.BF{1}.u/Results{2}.Data.BF{1}.norm.u,'LineWidth',2,'LineStyle',':','color',Color_map{2});
plot(DNS{1}.y(2,:,2)/delta_h,DNS{1}.u_avg_XZ./DNS{1}.u_b ,'LineWidth',2,'LineStyle','--','color',Color_map{3})
xlim([0 2])
ylim([0 1.5])
yticks([0 0.5 1.0 1.5])

xlabel('${y/\delta}$','interpreter','latex')
ylabel('$\langle u / U_r \rangle$','interpreter','latex')
pbaspect([1 1.0 1])

clengend = {'$HP$', '$IG$','$DNS$'};
legend(clengend,'interpreter','latex','fontsize',16,'Location','south',Box='off')

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)

exportgraphics(f,strcat('Figures/','Poiseuille_vs_DNS_streamwise','.jpeg'),'Resolution',300)



% T plot
f = figure
hold on; grid on; box on;

% Temperature
plot(Results{1}.Data.y,Results{1}.Data.BF{1}.T/Results{1}.Data.BF{1}.norm.T,'LineWidth',2,'LineStyle','-','color',Color_map{1});
plot(Results{2}.Data.y,Results{2}.Data.BF{1}.T/Results{2}.Data.BF{1}.norm.T,'LineWidth',2,'LineStyle',':','color',Color_map{2});
plot(DNS{1}.y(2,2:end-1,2)/delta_h,DNS{1}.T_avg_XZ(2:end-1)./DNS{1}.T_b ,'LineWidth',2,'LineStyle','--','color',Color_map{3})
xlim([0 2])
ylim([0.75 1.5])
yticks([0.75 1.0 1.25 1.5])

xlabel('${y/\delta}$','interpreter','latex')
ylabel('$\langle T / T_b \rangle$','interpreter','latex')
pbaspect([1 1.0 1])

clengend = {'$HP$', '$IG$','$DNS$'};
legend(clengend,'interpreter','latex','fontsize',16,'Location','northwest',Box='off')

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)

exportgraphics(f,strcat('Figures/','Poiseuille_vs_DNS_T','.jpeg'),'Resolution',300)


% C_p plot
f = figure
hold on; grid on; box on;
plot(Results{1}.Data.y,Results{1}.Data.BF{1}.c_p/Results{1}.Data.BF{1}.norm.c_p,'LineWidth',2,'LineStyle','-','color',Color_map{1});
plot(Results{2}.Data.y,Results{2}.Data.BF{1}.c_p/Results{2}.Data.BF{1}.norm.c_p,'LineWidth',2,'LineStyle',':','color',Color_map{2});
plot(DNS{1}.y(2,2:end-1,2)/delta_h,DNS{1}.c_p_avg_XZ(2:end-1)./DNS{1}.c_p_b ,'LineWidth',2,'LineStyle','--','color',Color_map{3})
xlim([0 2])
ylim([0.5 2.0])
% yticks([0 0.4 0.8 1.2 1.6])

xlabel('${y/\delta}$','interpreter','latex')
ylabel('$\langle c_p / {c_p}_b \rangle$','interpreter','latex')
pbaspect([1 1.0 1])

clengend = {'$HP$', '$IG$','$DNS$'};
legend(clengend,'interpreter','latex','fontsize',16,'Location','northwest',Box='off')

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)

exportgraphics(f,strcat('Figures/','Poiseuille_vs_DNS_c_p','.jpeg'),'Resolution',300)
