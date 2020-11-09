function out_RyR = get_h_RyR_reaction(h_RyR, cav,param)

h_infty = (param.K_hRyR^2)./(cav.^2+param.K_hRyR.^2);
out_RyR = (h_infty - h_RyR)/param.tau;


end

