function reactions = getborderreactions(c,ip,...
        ce,h,param, w_IPR,rate)



    
    % Calcium
    gamma = param.gamma*rate;
    % IPR
    k_beta = param.k_beta;
    K_p = param.K_p;
    K_c = param.K_c ;
    K_h = param.K_h;
    
    phi_p=ip.^2./(K_p^2+ip.^2);
    phi_c=c.^4./(c.^4+K_c^4);
    h_alpha = K_h^4./(c.^4+K_h^4);

    beta=phi_c.*phi_p.*h;
    alpha = (1-phi_p).*(1-phi_c.*h_alpha);
    P0 = beta./(beta+k_beta*(beta+alpha));

    JIPR=w_IPR.*P0;
    
    reactions(:,1) = (JIPR).*(ce-c);


    
end

