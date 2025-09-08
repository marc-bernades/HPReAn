function plot_TripletsEvol(File_sweep)

for ff = 1:length(File_sweep)
    % Unpack structure and load variables
    Results{ff} = load(strcat('Results/', File_sweep{ff}, '.mat'));
end

% Wavenumbers plots
for ff = 1:2

    % Compute Sigma
    ii = 1; jj = 1; cc = 1; % First and only position Re, Br and c
    for bb = 1:length(Results{ff}.Data.beta)
        for aa = 1:length(Results{ff}.Data.alpha)
            Results{ff}.Sigma_1(aa,bb) = Results{ff}.Data.S{aa,bb,ii,jj,cc}(1,1);
            Results{ff}.Sigma_2(aa,bb) = Results{ff}.Data.S{aa,bb,ii,jj,cc}(2,2);
            Results{ff}.Sigma_3(aa,bb) = Results{ff}.Data.S{aa,bb,ii,jj,cc}(3,3);

        end
    end

    % Streamwise dependency
    beta_target = 2.0;
    [~,bb] = min(abs(Results{ff}.Data.beta-beta_target));
    beta_target = Results{ff}.Data.beta(bb);
    Sigma_1_alpha{ff} = Results{ff}.Sigma_1(:,bb);
    Sigma_2_alpha{ff} = Results{ff}.Sigma_2(:,bb);
    Sigma_3_alpha{ff} = Results{ff}.Sigma_3(:,bb);

    % Spanwise dependency
    alpha_target = 0.0;
    [~,aa] = min(abs(Results{ff}.Data.alpha-alpha_target));
    alpha_target = Results{ff}.Data.alpha(aa);
    Sigma_1_beta{ff} = Results{ff}.Sigma_1(aa,:);
    Sigma_2_beta{ff} = Results{ff}.Sigma_2(aa,:);
    Sigma_3_beta{ff} = Results{ff}.Sigma_3(aa,:);

end

f = figure;
hold on; grid on; box on;
semilogx(Results{1}.Data.alpha, log10(Sigma_1_alpha{1}),'LineWidth',2,'LineStyle','-','color',[0.8500 0.3250 0.0980]);
semilogx(Results{1}.Data.alpha, log10(Sigma_2_alpha{1}),'LineWidth',2,'LineStyle','--','color',[0.8500 0.3250 0.0980]);
semilogx(Results{1}.Data.alpha, log10(Sigma_3_alpha{1}),'LineWidth',2,'LineStyle','-.','color',[0.8500 0.3250 0.0980]);
semilogx(Results{2}.Data.alpha, log10(Sigma_1_alpha{2}),'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
semilogx(Results{2}.Data.alpha, log10(Sigma_2_alpha{2}),'LineWidth',2,'LineStyle','--','color',[0 0.4470 0.7410])
semilogx(Results{2}.Data.alpha, log10(Sigma_3_alpha{2}),'LineWidth',2,'LineStyle','-.','color',[0 0.4470 0.7410])

xlabel('${k_x}$','interpreter','latex')
ylabel('${{log}_{10}(\sigma_{1})}$','interpreter','latex')

h(1) = plot( NaN,'-','color','k','LineWidth',2);
h(2) = plot( NaN,'--','color','k','LineWidth',2);
h(3) = plot(NaN,'-.','color','k','LineWidth',2);
idx_add = length(h);

setup_name = {'$HP$','$IG$'};
for kk = 1:2
    if kk == 1
        col_plot = [0.8500 0.3250 0.0980];
    else
        col_plot = [0 0.4470 0.7410];
    end

    h(kk+idx_add) = plot( NaN,'-','color',col_plot,'LineWidth',2);
    clegend{kk}   = strcat('$', setup_name{kk}, '$');
end
% clengend_final = {'$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|T^{\prime}|$',clegend{:}};
clengend_final = {'$\sigma_1$','$\sigma_2$','$\sigma_3$',clegend{:}};

legend(h,clengend_final,'interpreter','latex','NumColumns',2,'fontsize',16,'Location','northeast',Box='off')



% legend('$HP$','$IG$','interpreter','latex')
% legend('Location','northeast','box','off')

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)

set(gca,'XScale', 'log')


exportgraphics(f,strcat('Figures/','TripletEvol_HP_vs_IG_alpha_sigma1_3','.jpeg'),'Resolution',300)


f = figure;
hold on; grid on; box on;
semilogx(Results{1}.Data.beta, log10(Sigma_1_beta{1}),'LineWidth',2,'LineStyle','-','color',[0.8500 0.3250 0.0980]);
semilogx(Results{1}.Data.beta, log10(Sigma_2_beta{1}),'LineWidth',2,'LineStyle','--','color',[0.8500 0.3250 0.0980]);
semilogx(Results{1}.Data.beta, log10(Sigma_3_beta{1}),'LineWidth',2,'LineStyle','-.','color',[0.8500 0.3250 0.0980]);
semilogx(Results{2}.Data.beta, log10(Sigma_1_beta{2}), 'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
semilogx(Results{2}.Data.beta, log10(Sigma_2_beta{2}), 'LineWidth',2,'LineStyle','--','color',[0 0.4470 0.7410])
semilogx(Results{2}.Data.beta, log10(Sigma_3_beta{2}), 'LineWidth',2,'LineStyle','-.','color',[0 0.4470 0.7410])
xlabel('${k_z}$','interpreter','latex')
ylabel('${{log}_{10}(\sigma_{1})}$','interpreter','latex')

% legend('$HP$','$IG$','interpreter','latex')
% legend('Location','northwest','box','off')
h = [];
h(1) = plot( NaN,'-','color','k','LineWidth',2);
h(2) = plot( NaN,'--','color','k','LineWidth',2);
h(3) = plot(NaN,'-.','color','k','LineWidth',2);
idx_add = length(h);

setup_name = {'$HP$','$IG$'};
for kk = 1:2
    if kk == 1
        col_plot = [0.8500 0.3250 0.0980];
    else
        col_plot = [0 0.4470 0.7410];
    end

    h(kk+idx_add) = plot( NaN,'-','color',col_plot,'LineWidth',2);
    clegend{kk}   = strcat('$', setup_name{kk}, '$');
end
% clengend_final = {'$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|T^{\prime}|$',clegend{:}};
clengend_final = {'$\sigma_1$','$\sigma_2$','$\sigma_3$',clegend{:}};

legend(h,clengend_final,'interpreter','latex','NumColumns',2,'fontsize',16,'Location','southwest',Box='off')



% legend('$HP$','$IG$','interpreter','latex')
% legend('Location','northeast','box','off')

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)

set(gca,'XScale', 'log')

exportgraphics(f,strcat('Figures/','TripletEvol_HP_vs_IG_beta_sigma1_3','.jpeg'),'Resolution',300)


% Frequency plots
for ff = 3:4

    % Compute Sigma
    ii = 1; jj = 1; aa = 1; bb = 1; % First and only position Re, Br, alpha and beta
    for cc = 1:length(Results{ff}.Data.c)
        Results{ff}.Sigma_1(cc) = Results{ff}.Data.S{aa,bb,ii,jj,cc}(1,1);
        Results{ff}.Sigma_2(cc) = Results{ff}.Data.S{aa,bb,ii,jj,cc}(2,2);
        Results{ff}.Sigma_3(cc) = Results{ff}.Data.S{aa,bb,ii,jj,cc}(3,3);
    end

    % Frequency dependency
    Sigma_1_omega{ff-2} = Results{ff}.Sigma_1;
    Sigma_2_omega{ff-2} = Results{ff}.Sigma_2;
    Sigma_3_omega{ff-2} = Results{ff}.Sigma_3;

    omega{ff-2} = Results{ff}.Data.c.*Results{ff}.Data.alpha;

end

% Frequency dependency > Equivalent to alpha

f = figure;
hold on; grid on; box on;
semilogx(omega{1}, log10(Sigma_1_omega{1}),'LineWidth',2,'LineStyle','-','color',[0.8500 0.3250 0.0980]);
semilogx(omega{1}, log10(Sigma_2_omega{1}),'LineWidth',2,'LineStyle','--','color',[0.8500 0.3250 0.0980]);
semilogx(omega{1}, log10(Sigma_3_omega{1}),'LineWidth',2,'LineStyle','-.','color',[0.8500 0.3250 0.0980]);
semilogx(omega{2}, log10(Sigma_1_omega{2}), 'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
semilogx(omega{2}, log10(Sigma_2_omega{2}), 'LineWidth',2,'LineStyle','--','color',[0 0.4470 0.7410])
semilogx(omega{2}, log10(Sigma_3_omega{2}), 'LineWidth',2,'LineStyle','-.','color',[0 0.4470 0.7410])
xlabel('$\omega$','interpreter','latex')
ylabel('${{log}_{10}(\sigma_{1})}$','interpreter','latex')

% legend('$HP$','$IG$','interpreter','latex')
% legend('Location','northeast','box','off')

h = [];
h(1) = plot( NaN,'-','color','k','LineWidth',2);
h(2) = plot( NaN,'--','color','k','LineWidth',2);
h(3) = plot(NaN,'-.','color','k','LineWidth',2);
idx_add = length(h);

setup_name = {'$HP$','$IG$'};
for kk = 1:2
    if kk == 1
        col_plot = [0.8500 0.3250 0.0980];
    else
        col_plot = [0 0.4470 0.7410];
    end

    h(kk+idx_add) = plot( NaN,'-','color',col_plot,'LineWidth',2);
    clegend{kk}   = strcat('$', setup_name{kk}, '$');
end
% clengend_final = {'$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$', '$|T^{\prime}|$',clegend{:}};
clengend_final = {'$\sigma_1$','$\sigma_2$','$\sigma_3$',clegend{:}};

legend(h,clengend_final,'interpreter','latex','NumColumns',2,'fontsize',16,'Location','southwest',Box='off')



% legend('$HP$','$IG$','interpreter','latex')
% legend('Location','northeast','box','off')

set(gca,'linewidth',1.5)
set(gca,'fontsize',16)

set(gca,'XScale', 'log')

exportgraphics(f,strcat('Figures/','TripletEvol_HP_vs_IG_omega_sigma1_3','.jpeg'),'Resolution',300)



end