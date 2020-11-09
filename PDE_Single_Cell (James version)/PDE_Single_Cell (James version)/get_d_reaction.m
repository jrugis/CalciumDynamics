function reactions= get_d_reaction(c,h,param)
%% Description

%% Input

%% Output

%% 
%     K3K=0.4;
%     Vdeg=0.1; gamma=0.185; Vs=0.25; Ks=0.1;
%     e=1;  kleak=0.00148; kRyR=0.01;
%     Ka=0.3; Ki=0.06; Kp=0.5;
%     tau=0.5;
%     K_RyR=0.42;

K_h=param.K_h;
tau_max = param.tau_max;
K_tau = param.K_tau;
h_inf=K_h^4./(K_h^4+c.^4);
tau =tau_max*K_tau^4./(K_tau^4+c.^4);

reactions = (h_inf-h)./tau;

end

