clear;close all;clc;format longg;load('Par.mat')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options=odeset('RelTol',1e-12,'AbsTol',1e-12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ca= fca(par,0.0751,0.0751);
Ca= fca(par,0.1739,0.21739);
f_secretion = @(t,x) cluster(t,x,Ca,par);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t,U] = ode15s(f_secretion,[0 5e2],par.IC_new, options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w = {U(:,1),U(:,9),U(:,17),U(:,25),U(:,33),U(:,41),U(:,49)};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
for i = 1:length(w)
    plot(t,w{i},'LineWidth',2)
    hold on
end
ax = gca;
ax.LineWidth = 2.1;
ax.FontSize = 20;
xlabel('Time (sec)')
ylabel('Volume (\mum^3)')
box off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run FFR

figure(2)
plot(t,Q33,'LineWidth',2)
hold on
ax = gca;
ax.LineWidth = 2.1;
ax.FontSize = 20;
xlabel('Time (sec)')
ylabel('Cluster Flow Rate (\mum^3/sec)')
box off


clear