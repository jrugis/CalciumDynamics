
close all
clear all
clc

param.VPLC = 0.01;
param.VPLCkeep = param.VPLC;  % Needed for setting VPLC = 0 during the first 200 s
% Parameters of diffusion
Dc=5.0;     % Calcium
Dp=283;%300.0;   % Ip3
param.Dp = Dp;
De=5.0;     % Calcium in the ER

% Parameters of transmission. Because this is for a single cell, these
% parameter are not used
%Fc = 0;        % intercellular transmission parameter Calcium
%Fip = 0;       % intercellular transmission parameter Ip3

% Parameters for time scheme
delt=0.1;  % Time step
param.initialperiod = 200;   % time you let the cell run before stimulating. Should be about 200 for steady state.
param.stimperiod = 60;      % length of the stimulation period
tend = param.initialperiod + param.stimperiod;   % Total time
t_keep = tend;

cells_to_simulate = 1:1:2;
param.cells_to_simulate = cells_to_simulate;
n_cells = max(size(cells_to_simulate));

%V_RyR = 0.15;
param.V_RyR= 0;   % This param should be set to zero in the new model
%K_RyR = [0.15,0.3,0.42];
param.K_RyR = 0.2;  % 0.16 for more sensitive RyR
param.K_RyR2 = 17;
param.d_RyR = 1;
param.n_RyR = 4;
param.m_RyR = 0;

% IPR
param.k_beta =0.4;
param.K_p = 0.2;
param.K_c = 0.2;
param.K_h = 0.08;
param.kIPR = 45;
param.kIPR_min = 0;

% Serca
param.V_p = 0.9;
param.k_p = 0.2;
param.K_bar = 0.00001957;

% Calcium
param.gamma=0.185;

% IP3
param.V_3K=0.05;
param.V_5K=0.05;
param.K_PLC = 0.07;
param.K3K=0.4;

% h
param.tau_max = 100;
param.K_tau = 0.1;

% h RyR
param.K_hRyR = 0.15;
param.tau = 1;                    

param.IPRdn = 0.8;      % Width of band of apical IPR.
param.IPRdf = 0.11;     % This doesn't seem to be used anywhere. Not sure what it does.
param.PLCds = 0.8;
param.PLCdl = 1.2;

% Ca2+-activated-K/Cl channels
param.KCaCC =0.26;     % uM
param.KCaKC =0.35;      % uM
param.eta1  =4.49;      % Hill Coefficient
param.eta2  =2;       % Hill Coefficient

% NaK
param.r = ( 1.305e6 ) *10^(-9);       % mM^-3 s^-1

param.alpha1 = ( 0.641 ) * 10^(-3) ;  % mM^-1

% Thermodynamical Constants
param.R     =8.3144621;    % J mol^-1 K^-1
param.T     =310;          % K
param.F     =96485.3365;   % C mol^-1
param.RTF   =1e3*param.R*param.T/param.F;

param.s = 1e-3;

% Steady state Values for pH - Pena-Muntzenmayer et al. (2015)
param.pHl = 6.81;
param.pHi = 6.91;
param.pHe = 7.41;

% Bicarbonate Buffer
param.kn = 2.64e4 * 0.012;
param.kp = 11 * 0.012;


% Na-K-2Cl
param.a1 = 157.55;                  %       s^-1
param.a2 = (2.0096*1e7)*(10^-12);   % mM^-4 s^-1
param.a3 = 1.0306;                        % s^-1
param.a4 = 1.3852*1e6*(10^-12);     % mM^-4 s^-1

% Steady state Concentrations
% Interstitum
param.Cle   =102.6;          %mM
param.Nae   =140.2;          %mM
param.Ke    =5.3;             %mM
param.He    =1e3*10^(-param.pHe);     %mM
param.HCO3e =40;
param.CO2e  =1.6;

% Anion Exchangers Falkenberg et al. (2010)
param.KCl = 5.6;
param.KB = 100;
param.KNa = 15;
param.KH = 4.5e-4;

% Membrane permeability to water
param.B1 = 0.0443354480372293 * 100;
param.B2 = 0.0194580852061638 * 100;
param.B3 = 0.000593663800587253 * 100;

c0_init = 0.0717;
h_init = param.K_h^4./(param.K_h^4+c0_init.^4);
param.ct = 3;

c0 = [c0_init,c0_init*ones(1,n_cells-1)]; % Initial condition Calcium in each cell
ip_init = 0;%.5*(2e-4/0.1)*((param.K3K^2+c0_init^2)/(c0_init^2));
ip0 = [ip_init,ip_init*ones(1,n_cells-1)]; % Initial condition Ip3 in each cell
ce0 = ((param.ct-c0)/param.gamma);            % Initial condition ER Calcium in each cell
h0 = 0*h_init*ones(1,n_cells);
hRyR0 = 0*c0.^2./(c0.^2+param.K_hRyR^2);

Nal0=118.7;             % mM
Kl0=5.6;                % mM
Cll0=Nal0+Kl0;          % mM
Hl = 1e3*10^(-param.pHl);   %mM
HCO3l = Kl0 + Nal0 - Cll0 + Hl;
CO2l = 11.6;

Na0=25;                 % mM
K0=120;                 % mM
Cl0=50;                 % mM
H0 = 1e3*10^(-param.pHi);     % mM
HCO30 = 10;             % mM
param.CO20 = ( 0.197e4 * ( CO2l + param.CO2e ) - param.kn * HCO30 * H0 ) /...
    ( 2 * 0.197e4 - param.kp );
param.Ul = ( param.B2 / param.B1 ) * ( 2 * ( Na0 + K0 + H0 ) + param.CO20 - ...
    ( param.Nae + param.Ke + param.Cle + param.HCO3e ) ) - ...
    ( 2 * ( Nal0 + Kl0 - Na0 - K0 - H0 ) - param.CO20 );

% Membrane potential SS
Va0=-50.24;         % mV
Vb0=-62.8;          % mV
Vt0=Va0-Vb0;        % mV

VtNa0 = param.RTF * log( Nal0 / param.Nae);
VtK0 = param.RTF * log( Kl0 / param.Ke );
VCl0 = param.RTF * log( Cll0 / Cl0);
VK0 = param.RTF * log( param.Ke / K0 );
CO20 = ( 0.197e4 * ( CO2l + param.CO2e ) -...
    param.kn * HCO30 * H0 ) / ( 2 * 0.197e4 - param.kp );
Qa0 = 16.540088976029931;
Qt0 = 0.726112417038916;
Qtot0 = 17.266201393068847;
Qb0 = 16.540088976029921;

file_mesh = {'out_N4_p3-p2-p4-1tet.msh','out_N4_p3-p2-p4-2tet.msh','out_N4_p3-p2-p4-3tet.msh',...
    'out_N4_p3-p2-p4-4tet.msh','out_N4_p3-p2-p4-5tet.msh',...
    'out_N4_p3-p2-p4-6tet.msh','out_N4_p3-p2-p4-7tet.msh'};

show = 0;

mod_basal = 1;

save('single_cell_params')







