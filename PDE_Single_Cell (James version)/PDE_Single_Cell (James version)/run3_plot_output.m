 clc; 
 close all; format long; 
 clear all;
 
 load('mod_basal1data_smoothed_mesh.mat')
 load('single_cell_params.mat')
 load('single_cell_conductance_file.mat')

     %load('c4_old_RyR_VPLC0.01.mat')
     %load('c4_old_noRyR_VPLC0.01.mat') 
     %load('c4_new_VPLC0.01.mat')
load('single_cell_output.mat')
 
chop = 15;               % chop off the last few points which are always weird, for some unknown reason.
ion_tot{1} =  ion_tot{1}(:,1:end-chop);
reg_time = reg_time(1:end-chop);
c_tot{1} = c_tot{1}(:,1:end-chop);

for i = 1:1
    Nal     = ion_tot{i}(1,:);
    Kl 		= ion_tot{i}(2,:);  
    Cl 		= ion_tot{i}(3,:);
    vol     = ion_tot{i}(4,:);
    Na 		= ion_tot{i}(5,:);
    K 		= ion_tot{i}(6,:);
    C 		= ion_tot{i}(7,:);
    H		= ion_tot{i}(9,:);
    Va      = ion_tot{i}(10,:);
    Vb      = ion_tot{i}(11,:);
    
    Qa = param.B1 * ( 2 * ( Nal + Kl - Na - K - H ) - param.CO20 + param.Ul );     % micro-metres^3.s^-1
    Qb = param.B2 * ( 2 * ( Na + K + H ) + param.CO20 - ...
                      ( param.Nae + param.Ke + param.Cle + param.HCO3e ) );
    Qt = param.B3 * ( 2 * ( Nal + Kl ) + param.Ul - ....
                      ( param.Nae + param.Ke + param.Cle + param.HCO3e ) ); % micro-metres^3.s^-1
    Qtot{i}=(Qa+Qt);                                     % micro-metres^3.s^-1
end
QFFR = Qtot{1};

ncell =2;  % Which cell was the simulation done on?

dum  = randi(size(p{ncell,1},1),20,1);  % select 20 random grid points for plotting calcium
dum = find(dist_ap_p{ncell,1}<0.2);     % select the apical grid points only for plotting calcium



%% Panel of fluid flow and ion plots. And calcium.
h2=figure(2);

subplot(2,3,6)
plot(reg_time,c_tot{1}(dum,:))
xlim([200 259])

subplot(2,3,1)
plot(reg_time,QFFR,'b','LineWidth',1.5)
box off
xlabel('Time')
ylabel('FFR')
xlim([200 259])

subplot(2,3,2)
plot(reg_time,vol,'LineWidth',1.5)
legend('cell volume')
xlim([200 259])

subplot(2,3,3)
plot(reg_time,Nal,reg_time,Kl,reg_time,Cl,'LineWidth',1.5)
legend('Na_l','K_l','Cl_l')
xlim([200 259])

subplot(2,3,4)
plot(reg_time,Na,reg_time,K,reg_time,C,'LineWidth',1.5)
legend('Na','K','Cl')
xlim([200 259])

subplot(2,3,5)
plot(reg_time,Va,reg_time,Vb,'LineWidth',1.5)
legend('Va','Vb')

set(h2,'Position',[1400 600 1350 700]);

%% Scatter plot of the points where you plot the Ca outputs
h4 = figure(4);
X=p{ncell,1}(dum,1);
Y=p{ncell,1}(dum,2);
Z=p{ncell,1}(dum,3);
scatter3(X,Y,Z,30,'b')
set(h4,'Position',[1900 100 500 400]);


%% Some test stuff. Plots the apical tets and triangles.
 h5=figure(5);
 tetramesh(tets{ncell}(find(dist_ap{ncell}<0.6),:),p{ncell});
 set(h5,'Position',[1400 100 500 400]);
 
 h6 = figure(6);
 dum=mean(dist_ap_p{ncell}(triangles{ncell}),2);
 trisurf(triangles{ncell}(find(dum<0.6),:),p{ncell}(:,1),p{ncell}(:,2),p{ncell}(:,3))
 set(h6,'Position',[1000 100 500 400]);

%% output files for Igor. The renaming is for convenience inside Igor. You have to comment out the bits
%     you don't want as I'm too lazy to write it generally. So sue me.

% old_noRyR_cal = c_tot{1}(dum,:);
% old_noRyR_QFFR = QFFR;
% old_noRyR_vol = vol;
% old_noRyR_time = reg_time;
% save('forigorplots_old_noRyR.mat','old_noRyR_time','old_noRyR_QFFR','old_noRyR_cal','old_noRyR_vol')

% old_RyR_cal = c_tot{1}(dum,:);
% old_RyR_QFFR = QFFR;
% old_RyR_vol = vol;
% old_RyR_time = reg_time;
% save('forigorplots_old_RyR.mat','old_RyR_time','old_RyR_QFFR','old_RyR_cal','old_RyR_vol')

% new_cal = c_tot{1}(dum,:);
% new_QFFR = QFFR;
% new_vol = vol;
% new_time = reg_time;
% save('forigorplots_new.mat','new_time','new_QFFR','new_cal','new_vol')
     