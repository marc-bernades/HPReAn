function plot_ResponseForcing(y,delta,q_vars,U,S,V,aa,bb,ii,jj,oo, idx, alpha, beta, c, name_file_load)


for n_var = 1:q_vars
    Q{n_var}     = U{aa,bb,ii,jj,oo}(n_var:q_vars:end,idx); % Take all columns for each perturbation row
    F{n_var}     = V{aa,bb,ii,jj,oo}(n_var:q_vars:end,idx); % Take all columns for each perturbation row
end

Sigma(1,oo) = S{aa,bb,ii,jj,oo}(1,1);
Sigma(2,oo) = S{aa,bb,ii,jj,oo}(2,2);
Sigma(3,oo) = S{aa,bb,ii,jj,oo}(3,3);

disp("Sigma 1 = " + num2str(Sigma(1,oo)))
disp("Sigma 2 = " + num2str(Sigma(2,oo)))
disp("Sigma 3 = " + num2str(Sigma(3,oo)))

disp("")


%% Response
Norm = max(abs(Q{1})); % 1;
f = figure;
hold on; grid on; box on;
plot(y/delta,abs(Q{1}(:,idx))/Norm,'LineWidth',2,'LineStyle',':','color',[0.4660 0.6740 0.1880]);
plot(y/delta,abs(Q{2}(:,idx))/Norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
plot(y/delta,abs(Q{3}(:,idx))/Norm,'LineWidth',2,'LineStyle','--','color',[0.8500 0.3250 0.0980])
plot(y/delta,abs(Q{4}(:,idx))/Norm,'LineWidth',2,'LineStyle','-.','color','black')
plot(y/delta,abs(Q{5}(:,idx))/Norm,'LineWidth',2,'LineStyle','-.','color',[0.4940 0.1840 0.5560])
xlabel('${y/\delta}$','interpreter','latex')
ylabel('$|q^{\prime}|$','interpreter','latex')

legend('$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$','$|w^{\prime}|$', '$|T^{\prime}|$','interpreter','latex')
legend('Location','northeast','box','off')

set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

% Inset
% axes('position',[.25 .56 .35 .3])
% 
% hold on; grid off; box on;
% plot(y/delta,abs(Q{1}(:,idx))/Norm,'LineWidth',2,'LineStyle',':','color',[0.4660 0.6740 0.1880]);
% plot(y/delta,abs(Q{2}(:,idx))/Norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
% plot(y/delta,abs(Q{3}(:,idx))/Norm,'LineWidth',2,'LineStyle','--','color',[0.8500 0.3250 0.0980])
% plot(y/delta,abs(Q{4}(:,idx))/Norm,'LineWidth',2,'LineStyle','-.','color','black')
% plot(y/delta,abs(Q{5}(:,idx))/Norm,'LineWidth',2,'LineStyle','-.','color',[0.4940 0.1840 0.5560])
% 
% ylim([0 0.01])
% 
% set(gca,'linewidth',1.5)
% set(gca,'fontsize',10)

exportgraphics(f,strcat('Figures/Response_',name_file_load,'_c_',num2str(c*alpha(aa)),...
    '_alpha_',num2str(alpha(aa)),'_beta_',num2str(beta(bb)),'.jpeg'),'Resolution',300)


%% Forcing
Norm = max(abs(F{3})); % 1;
f = figure;
hold on; grid on; box on;
plot(y(3:2:end-2)/delta,abs(F{1}((3:2:end-2),idx))/Norm,'LineWidth',2,'LineStyle',':','color',[0.4660 0.6740 0.1880]);
plot(y/delta,abs(F{2}(:,idx))/Norm,'LineWidth',2,'LineStyle','-','color',[0 0.4470 0.7410])
plot(y/delta,abs(F{3}(:,idx))/Norm,'LineWidth',2,'LineStyle','--','color',[0.8500 0.3250 0.0980])
plot(y/delta,abs(F{4}(:,idx))/Norm,'LineWidth',2,'LineStyle','-.','color','black')
plot(y/delta,abs(F{5}(:,idx))/Norm,'LineWidth',2,'LineStyle','-.','color',[0.4940 0.1840 0.5560])
xlabel('${y/\delta}$','interpreter','latex')
ylabel('$|f^{\prime}|$','interpreter','latex')
%
legend('$|\rho^{\prime}|$','$|u^{\prime}|$','$|v^{\prime}|$','$|w^{\prime}|$','$|T^{\prime}|$','interpreter','latex')
legend('Location','northeast','box','off')

set(gca,'linewidth',1.5)
set(gca,'fontsize',14)

exportgraphics(f,strcat('Figures/Forcing_',name_file_load,'_c_',num2str(c*alpha(aa)),...
    '_alpha_',num2str(alpha(aa)),'_beta_',num2str(beta(bb)),'.jpeg'),'Resolution',300)


end





