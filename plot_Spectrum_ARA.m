function plot_Spectrum_ARA(y, delta, N, alpha, N_target,n_LST_Sweep, Val_eigen, name_file_load, x_lim, y_lim)

realc = real(Val_eigen);
imagc = imag(Val_eigen);
disp(['Most unstable mode = ', num2str(max(imagc)), ' total of = ', num2str(sum(imagc>0))])

figure
hold on; grid on; box on
plot(realc,imagc,'o','color',[0 0.4470 0.7410],'MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410])
% plot(imagc,realc,'o','color',[0 0.4470 0.7410],'MarkerSize',5,'MarkerFaceColor',[0 0.4470 0.7410])

% xlim(x_lim)
% ylim(y_lim)

xlabel('${\omega}_r$','interpreter','latex')
ylabel('${\omega}_i$','interpreter','latex')

set(gca,'fontsize',16)
set(gca,'linewidth',1.5)




% saveas(gca,strcat('Figures/', name_file_load,'_Spectrum_PrEc_',num2str(N_target),'_alpha_',num2str(alpha), '_Re_',num2str(n_LST_Sweep),'_N_', num2str(N)),'epsc')
% saveas(gca,strcat('Figures/', name_file_load,'_Spectrum_PrEc_',num2str(N_target),'_alpha_',num2str(alpha), '_Re_',num2str(n_LST_Sweep),'_N_', num2str(N)),'png')

end