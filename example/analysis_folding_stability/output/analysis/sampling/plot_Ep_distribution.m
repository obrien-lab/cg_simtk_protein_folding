clear all
n_skip = 20000;
T_windows = [306.0 316.0 326.0 333.0 339.0 343.0 345.0 346.0 347.0 349.0 353.0 359.0 367.0 377.0 387.0];
legend_str = {};
for i=1:length(T_windows)
    legend_str{i} = ['T = ' num2str(T_windows(i)) ' K'];
end

figure(1)
set(gcf,'Units','centimeters','Position',[8 1.5 12 10],'paperpositionmode','auto');
C = colormap(jet(length(T_windows)));
hold on
for i=1:length(T_windows)
    D = [];
    D = load(['../../aa' num2str(i) '/ene_1.log']);
    [N, edges] = histcounts(D(n_skip+1:end),50,'Normalization','probability');
    x = (edges(1:end-1)+edges(2:end))/2;
    patch(x,N,C(i,:),'FaceAlpha',0.4,'EdgeColor','k','LineWidth',1.0)
end
set(gca, 'fontsize',11,'fontweight','normal','LineWidth',1.5,'fontname','Nimbus Roman No9 L')
box on
grid on
xlabel('Potential Energy (kcal/mol)','fontsize',13,'color','k','Interpreter','latex')
ylabel('Probability','fontsize',13,'color','k','Interpreter','latex')
h = legend('String',legend_str,'Location','best','fontsize',6,'box','off','Interpreter','latex');
saveas(gcf, 'Ep_distribution.svg');
quit