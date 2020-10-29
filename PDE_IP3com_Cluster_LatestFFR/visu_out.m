function  visu_out(file,folder,n_cells)
if folder
load([num2str(folder),'/',file])
load([num2str(folder),'/','sol',file])
%load('conductancesFc0Fip0t1nc148c4e7af77f29f59f4cacdebb8053cf4.mat')
else
    load(file)
    load(['sol',file])   
    %load('conductancesFc0Fip0t1nc148c4e7af77f29f59f4cacdebb8053cf4.mat')
end
load('mod_basal1data_smoothed_mesh.mat')
%% Treat mesh to choose on which cell to act
p = select_cell(p,cells_to_simulate,n_cells);
tets = select_cell(tets,cells_to_simulate,n_cells);
volume = select_cell(volume,cells_to_simulate,n_cells);
tets_volume = select_cell(tets_volume,cells_to_simulate,n_cells);
membranes = select_cell(membranes,cells_to_simulate,n_cells);
triangles = select_cell(triangles,cells_to_simulate,n_cells);
bndry = select_cell(bndry,cells_to_simulate,n_cells);
dist_ap = select_cell(dist_ap,cells_to_simulate,n_cells);
dist_ap_p = select_cell(dist_ap_p,cells_to_simulate,n_cells);
dist_ba = select_cell(dist_ba,cells_to_simulate,n_cells);
%%
% per = [103,93,105,116,46,108,103];
per = [108,84,108,113,41,109];
time_t = (reg_time>=max(reg_time-600))&(reg_time<max(reg_time-1));
reg_time = reg_time(time_t);
for ncell =1:n_cells
    disp(num2str(cells_to_simulate(ncell)))
    tic


c_tot{ncell} = c_tot{ncell}(:,time_t);
ip_tot{ncell} = ip_tot{ncell}(:,time_t);
ion_tot= ion_tot(:,time_t);

ind_x = repmat((1:size(tets{ncell},1))',1,4);
ind_y = tets{ncell};
ind_x = ind_x(:);
ind_y = ind_y(:);
val = (1/4)*ones(4*size(tets{ncell},1),1);
m_tets = sparse(ind_x,ind_y,val,size(tets{ncell},1),size(p{ncell},1));

c_tets = m_tets*c_tot{ncell};
ip_tets = m_tets*ip_tot{ncell};
d = [0,0.8];
leg = cell(2,1);
% for i = 1:size(d,2)
%    leg{i}= [num2str(d(i)), '\mu m' ];
% end
leg{1} = 'Apical';
leg{2} = 'Basal';
c_vol = zeros(size(d,2),size(c_tot{ncell},2));

w_basal = (dist_ba{ncell}<0.8);
w_apical = (dist_ap{ncell}<0.8);

v_basal_tets = tets_volume{ncell}(w_basal);
v_apical_tets = tets_volume{ncell}(w_apical);
v_basal = sum(v_basal_tets);
v_apical = sum(v_apical_tets);

c_basal_tets = c_tets(w_basal,:);
c_apical_tets = c_tets(w_apical,:);
ip_basal_tets = ip_tets(w_basal,:);
ip_apical_tets = ip_tets(w_apical,:);

c_basal = sum(v_basal_tets'.*c_basal_tets)/v_basal;
c_vol = sum(v_apical_tets'.*c_apical_tets)/v_apical;
ip_basal = sum(v_basal_tets'.*ip_basal_tets)/v_basal;
ip_vol = sum(v_apical_tets'.*ip_apical_tets)/v_apical;

ip_total = sum(tets_volume{ncell}'.*ip_tets)/sum(tets_volume{ncell});

c_tot{ncell} = 0;

ip_tets = 0;
disp('Treated')
toc



%% Calcium in the cell
h = figure(1+(ncell-1)*4);
plot(reg_time,c_vol(1,:),'linewidth',5,'Color',[137,192,113]/255)
hold on
plot(reg_time,c_basal,'linewidth',5,'Color',[211,99,57]/255)
xlabel('Time (s)')
ylabel('[Ca^{2+}] \muM')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5625 1]);
legend(leg)
ax=gca;
set(ax,'Linewidth',10)
ax.FontSize=100;
box off
axis([0,300, 0.05, 0.4])
xticks([0 100 200 300])
yticks([0.1 0.3 0.5])

%% IP_3 in the cell
h = figure(2+(ncell-1)*4);
plot(reg_time,ip_vol(1,:),'linewidth',5,'Color',[137,192,113]/255)
hold on
plot(reg_time,ip_basal,'linewidth',5,'Color',[211,99,57]/255)
xlabel('Time (s)')
ylabel('[IP_3] \muM')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5625 1]);
legend(leg)
ax=gca;
set(ax,'Linewidth',10)
ax.FontSize=100;
box off
axis([0,300, 0.05, 0.4])
xticks([0 100 200 300])
yticks([0.1 0.3 0.5])

%%
h = figure(42);
hold on
plot(reg_time,ip_total,'linewidth',5)
xlabel('Time (s)')
ylabel('[IP_3] \muM')
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5625 1]);
legend(leg)
ax=gca;
set(ax,'Linewidth',10)
ax.FontSize=100;
box off
% axis([0,300, 0.05, 0.4])
% xticks([0 100 200 300])
% yticks([0.1 0.3 0.5])

%% Ions and flux


w = ion_tot(1,:);
Na = ion_tot(2,:);
K = ion_tot(3,:);
Cl = ion_tot(4,:);
HCO3 = ion_tot(5,:);
H = ion_tot(6,:);

% Qa = param.B1 * ( 2 * ( Nal + Kl - Na - K - H ) - param.CO20 + param.Ul );     % micro-metres^3.s^-1
% Qb = param.B2 * ( 2 * ( Na + K + H ) + param.CO20 - ...
% ( param.Nae + param.Ke + param.Cle + param.HCO3e ) );
% Qt = param.B3 * ( 2 * ( Nal + Kl ) + param.Ul - ...
% ( param.Nae + param.Ke + param.Cle + param.HCO3e ) ); % micro-metres^3.s^-1
% Qtot=(Qa+Qt);


figure(3+(ncell-1)*4),

subplot(4,3,1)
plot(reg_time,Na,'linewidth',3)
xlabel('Time (sec)')
ylabel('[Na^+]_i (mM)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on

subplot(4,3,2)
plot(reg_time,K,'linewidth',3)
xlabel('Time (sec)')
ylabel('[K^+]_i (mM)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on

subplot(4,3,3)
plot(reg_time,Cl,'linewidth',3)
xlabel('Time (sec)')
ylabel('[Cl^-]_i (mM)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on


subplot(4,3,6)
plot(reg_time,w,'linewidth',3)
xlabel('Time (sec)')
ylabel('Volume (\mum^3)')
ax=gca;
set(ax,'Linewidth',2)
ax.FontSize=20;
box off
hold on


% subplot(4,3,10)
% plot(reg_time,Qtot,'linewidth',3)
% xlabel('Time (sec)')
% ylabel('Water Flow (\mum^3/s)')
% ax=gca;
% set(ax,'Linewidth',2)
% ax.FontSize=20;
% box off
% hold on
% hold off
% 
% h = figure(4+(ncell-1)*4),
% plot(reg_time,Qtot,'linewidth',5)
% xlabel('Time (s)')
% ylabel('Water Flow (\mum^3/s)')
% hold on
% ax=gca;
% set(ax,'Linewidth',10)
% ax.FontSize=100;
% box off
% axis([0,300, 20, 160])
% xticks([0 100 200 300])
% yticks([20 100 160])
% set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.5625 1]);
% hold off
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(h,['FF_',num2str(cells_to_simulate(ncell)),'.pdf'],'-dpdf','-r0')
% dlmwrite('FF_PDE.txt',[reg_time;Qtot]',' ')
% [~,locs]=findpeaks(Qtot(end/2:end));
% diff(locs)
% av_per = mean(diff(locs)*0.1);
% av_FF = sum(Qtot(reg_time>(max(reg_time)-mean(diff(locs)*0.1))))*0.1/mean(diff(locs)*0.1);
% av_IP_ap = sum(ip_vol(1,reg_time>(max(reg_time)-mean(diff(locs)*0.1))))*0.1/mean(diff(locs)*0.1);
% av_IP_ba = sum(ip_basal(reg_time>(max(reg_time)-mean(diff(locs)*0.1))))*0.1/mean(diff(locs)*0.1);
% av_Ca_ap = sum(c_vol(1,reg_time>(max(reg_time)-mean(diff(locs)*0.1))))*0.1/mean(diff(locs)*0.1);
% max_Ca_ap = max(c_vol(1,reg_time>(max(reg_time)-mean(diff(locs)*0.1))));
% min_Ca_ap = min(c_vol(1,reg_time>(max(reg_time)-mean(diff(locs)*0.1))));
% av_Ca_ba = sum(c_basal(reg_time>(max(reg_time)-mean(diff(locs)*0.1))))*0.1/mean(diff(locs)*0.1);
% max_Ca_ba = max(c_basal(reg_time>(max(reg_time)-mean(diff(locs)*0.1))));
% min_Ca_ba = min(c_basal(reg_time>(max(reg_time)-mean(diff(locs)*0.1))));
% disp(['This simulation has an average period of ',num2str(av_per),'seconds,\n'])
% disp([' an average Fluid flow of ',num2str(av_FF),'\mu m^3 .s^-1,\n'])
% disp([' an average IP_3 concentration in the apical region of ',num2str(av_IP_ap),'\muM,'])
% disp([' an average IP_3 concentration in the basal region of ',num2str(av_IP_ba),'\muM,'])
% disp([' an average Ca^2+ concentration in the apical region of ',num2str(av_Ca_ap),'\muM,'])
% disp([' an Maximum Ca^2+ concentration in the apical region of ',num2str(max_Ca_ap),'\muM,'])
% disp([' an Minimum Ca^2+ concentration in the apical region of ',num2str(min_Ca_ap),'\muM,'])
% disp([' an average Ca^2+ concentration in the basal region of ',num2str(av_Ca_ba),'\muM,'])
% disp([' an Maximum Ca^2+ concentration in the basal region of ',num2str(max_Ca_ba),'muM,'])
% disp([' an Minimum Ca^2+ concentration in the basal region of ',num2str(min_Ca_ba),'muM.'])
%     

end
end

