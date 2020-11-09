function [out,out_h,out_RyR] = make_load(u,h,h_RyR, p, param,w_IPR, w_VPLC, w_RyR, vol_tets,surf_tri,...
    m_tets,m_triangles,rate)

    % global stiff mass
    np = size(p,1);


    %   First np components of u    -   c
    %   Second np components of u   -   ip3
    %   Third np components of u    -   ce
    %   Fourth np components of u    -   h

    c=u(1:np);
    ip=u(np+1:2*np);
    ce=(param.ct*rate-c)*(1/param.gamma);


    cav = m_tets*c;
    ipav = m_tets*ip;
    ceav = m_tets*ce;
    hRyRav = m_tets*h_RyR;
    
    cav_tri  = m_triangles(:,1:size(h,1))*c(1:size(h,1));
    ipav_tri = m_triangles(:,1:size(h,1))*ip(1:size(h,1));
    ceav_tri = m_triangles(:,1:size(h,1))*ce(1:size(h,1));
    hav_tri  = m_triangles(:,1:size(h,1))*h;
    
    reactions = getrhsreactions(cav,ipav,ceav,hRyRav, param, w_VPLC,...
        w_RyR,rate);

    reactions_border = getborderreactions(cav_tri,ipav_tri,...
        ceav_tri,hav_tri,param, w_IPR,rate);

    out_RyR = get_h_RyR_reaction(h_RyR, c,param);

    
    load_c = (m_tets')*(vol_tets.*reactions(:,1));
    load_ip = (m_tets')*(vol_tets.*reactions(:,2));

    load_c  = load_c + (m_triangles')*(surf_tri.*reactions_border(:,1));
    
    load_h = get_d_reaction(c(1:size(h,1)),h,param);


    out(1:np,1)             =    load_c;       % Diffusing variable
    out(np+1:2*np,1)        =    load_ip;
    out_h      =    load_h;


    % The variables with no diffusion should be done differently. No loop over
    % the triangles is needed, just a simpler loop over the nodes.


end

