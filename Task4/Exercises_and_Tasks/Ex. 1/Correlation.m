function [ g ] = Correlation( a,b )
g=figure;
scatter(a,b,'linewidth',2)
hold on
range = [5 10];
plot(range,range,'k','linewidth',2)
[~,a,b,~,regression_string,...
    correlation_string,N_string] = Goodness_Of_Fit(a,b,'linear');
plot(range,a*range+b,'r','linewidth',2)
grid on
xlabel('\itv_{LOS} \rm[m s^{-1}] (SpinnerLidar)','fontsize',14)
ylabel('\itv_{LOS} \rm[m s^{-1}] (Centroid)','fontsize',14)
legend('Scatter','y = x',regression_string,'location','northwest')
text(9,6,correlation_string,'fontsize',14)
text(9,5.5,N_string,'fontsize',14)
set(gca,'fontsize',14)
set(gca,'position',get(gca,'position')+0.02*ones(1,4));
axis equal
end

