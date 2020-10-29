function reactions = getrhsreactions(c,ip,ce,h,param,w_VPLC,w_RyR,rate)


    % RyR
    V_RyR = param.V_RyR;
    K_RyR = param.K_RyR;
    K_RyR2 = param.K_RyR2;
    n_RyR = param.n_RyR;
    m_RyR = param.m_RyR;
    
    % h_RyR
    K_h = param.K_h;


    % Serca
    V_p = param.V_p ;
    k_p = param.k_p ;
    K_bar = param.K_bar;

    % Calcium
    gamma = param.gamma*rate;




    % IP3
    V_3K=param.V_3K;
    V_5K=param.V_5K;
    K_PLC=param.K_PLC;
    K3K=param.K3K;

    
    J_SERCA = V_p*(c.^2-K_bar*ce.^2)./(c.^2+k_p^2);
    J_RYR = w_RyR.*V_RyR.*(c.^n_RyR./(c.^n_RyR+K_RyR^n_RyR)).*(ce.^m_RyR./(ce.^m_RyR+K_RyR2^m_RyR)).*h;

    % Now define the ODE RHS

    reactions(:,1) = (J_RYR).*(ce-c)-J_SERCA;
    reactions(:,2) = w_VPLC.*c.^2./(c.^2+K_PLC^2)-(V_5K+V_3K*c.^2./(c.^2+K3K^2)).*ip;


end