

close all
clear all
clc

%% Loading the parameter and mesh files
load('single_cell_params.mat')
load('mod_basal1data_smoothed_mesh')
load('matrixFc.mat')


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

%% Parameters

% Parameters for time scheme

numt=t_keep/delt;

% Sizes
n_tets = zeros(n_cells,1);
n_triangles = zeros(n_cells,1);
np = zeros(n_cells,1);
for i= 1:n_cells
    np(i) = size(p{i},1);
    n_tets(i) = size(tets{i},1);
    n_triangles(i) = size(triangles{i},1);
end

%% Matrices from tri to tets
% Compute the volume of each tetrahedron for each mesh
vol_tets = cell(n_cells,1);

for i = 1:n_cells
    vol_tets{i} = zeros(n_tets(i),1);
    for j = 1:n_tets(i)
        vert = p{i}(tets{i}(j,:),:);
        J = vert(2:4,:)-repmat(vert(1,:),3,1);
        vol_tets{i}(j) = abs(det(J))/6;
    end
end

% Matrix from p to tets
m_tets = cell(n_cells,1);
for i = 1:n_cells
    ind_x = repmat((1:n_tets(i))',1,4);
    ind_y = tets{i};
    ind_x = ind_x(:);
    ind_y = ind_y(:);
    val = (1/4)*ones(4*n_tets(i),1);
    m_tets{i} = sparse(ind_x,ind_y,val,n_tets(i),np(i));
end


% Compute the surface of each triangle at the surface
surf_tri = cell(n_cells,1);

for i = 1:n_cells
    surf_tri{i} = zeros(n_triangles(i),1);
    for j = 1:n_triangles(i)
        vert = p{i}(triangles{i}(j,:),:);
        J = vert(2:3,:)-repmat(vert(1,:),2,1);
        surf_tri{i}(j) = norm(cross(J(1,:),J(2,:)))/2;
    end
end

% Matrix from p to triangles
m_triangles = cell(n_cells,1);
for i = 1:n_cells
    ind_x = repmat((1:n_triangles(i))',1,3);
    ind_y = triangles{i};
    ind_x = ind_x(:);
    ind_y = ind_y(:);
    val = (1/3)*ones(3*n_triangles(i),1);
    m_triangles{i} = sparse(ind_x,ind_y,val,n_triangles(i),np(i));
end



%% Distance to consider the basal and apical region

p_border = cell(n_cells,1);
dist_ap_border =  cell(n_cells,1);
dist_ap_tets = cell(n_cells,1);
dist_ba_border = cell(n_cells,1);
for i = 1:n_cells
    p_border{i} = p{i}(unique(triangles{i}),:);
    dist_ap_border{i} = sum(dist_ap_p{i}(triangles{i}),2)/3;
    dist_ap_tets{i} = sum(dist_ap_p{i}(tets{i}),2)/4;
    dist_ba_p = m_tets{i}'*dist_ba{i};
    dist_ba_border{i} = sum(dist_ba_p(triangles{i}),2)/3;
end

w_IPR = cell(n_cells,1);
w_VPLC  = cell(n_cells,1);
w_basal = cell(n_cells,1);
w_RyR = cell(n_cells,1);
for i = 1:n_cells
    w_VPLC{i} = param.VPLC*(dist_ba{i}< param.PLCds).*...
        (dist_ap{i}> param.PLCdl);
    w_IPR{i} = param.kIPR*(dist_ap_border{i} < param.IPRdn);
    w_basal{i}= (dist_ap_border{i}>param.PLCdl).*ismember(triangles{i},bndry{i},'rows');
    w_RyR{i} = dist_ap_tets{i}.*(dist_ap_tets{i}<param.d_RyR)+...
        1.*(dist_ap_tets{i}>param.d_RyR);
end


%% Choosing the matrices
Amat = select_cell(Amat,cells_to_simulate,n_cells);
mass = select_cell(mass,cells_to_simulate,n_cells);
stiff = select_cell(stiff,cells_to_simulate,n_cells);
sMem = select_cell(sMem,cells_to_simulate,n_cells);
sMass = select_cell(sMass,cells_to_simulate,n_cells);


%% Initialisation


% Initialisation of the volume variables that will diffuse
u = cell(n_cells,1);
for i = 1:n_cells
    u{i} = zeros(2*np(i),numt);
    u{i}(1:np(i),1)=c0(i);    
    u{i}((np(i)+1):2*np(i),1)=ip0(i);  
end

% Initialization of the h_ryR
h_RyR = cell(n_cells,1);
for i = 1:n_cells
    h_RyR{i} = hRyR0(i)*ones(np(i),numt);
end


% Initialisation of the ions variables
ion = cell(n_cells,1);

for i = 1:n_cells
    ion{i} = zeros(11,numt);
    IC=[ Nal0 Kl0 Cll0 volume(i) Na0 K0 Cl0 HCO30 H0 Va0 Vb0];
    ion{i}(:,1) = IC;
end

% Initialisation of the h border variable
h_border = cell(n_cells,1);
for i = 1:n_cells
    h_border{i} = zeros(size(unique(triangles{i}(:)),1),numt);
    h_border{i}(:,1) =  h0(i);
end

% To store the results
c_tot = cell(n_cells,1); 
ip_tot = cell(n_cells,1);
ion_tot = cell(n_cells,1);
reg_time = (tend-t_keep):0.1:tend;



%% Parameters for each cell
for i = 1:n_cells
    param.Sa(i) = sum(((w_IPR{i}>0)*1).*surf_tri{i});
    param.St(i) = param.Sa(i) / 0.943551157250391;                            % micro-metres^2
    param.Sb(i) = sum(((w_basal{i}>0)*1).*surf_tri{i});
    param.wl(i) = 0.02*volume(i);
    
    param.aNkcc1(i) = ( 0.0063812 ) * param.Sb(i);
    
    % Tight Junction Na current GtNa
    vtNa = param.St(i) * ( Vt0 - VtNa0 ) / param.F;
    param.GtNa(i) = Qtot0 * Nal0 / vtNa;
    
    
    % Tight junction K current GtK
    vtK = param.St(i) * ( Vt0 - VtK0 ) / param.F; 
    param.GtK(i) = Qtot0 * Kl0 / vtK;
    
    
    % Apical Ca2+ activated Cl channels GCl
    PCl = 1 ./ ( 1 + ( param.KCaCC ./ c0_init ).^param.eta1 );
    vCl = PCl * ( Va0 + VCl0) / param.F;
    param.GCl(i) = - Qtot0 * Cll0 / vCl;
    
    
    % Bicarbonate Buffer GB
    vBB = volume(i) * ( param.kp * CO20 - param.kn * HCO30 * H0 );
    param.GB(i) = 4251.79 / vBB;
    
    % Sodium Proton Antiporter G1
    JBB0 = param.GB(i) * vBB;
    vNhe1 = param.Sb(i) * ( ( param.Nae / ( param.Nae + param.KNa ) ) *...
            ( H0 / ( param.KH + H0 ) ) -...
            ( Na0 / ( Na0 + param.KNa ) ) * ...
            ( param.He / ( param.KH + param.He ) ) );
    param.G1(i)= JBB0 / vNhe1;
    
    % Sodium Potassium ATPase aNaK
    JtNa0 = param.GtNa(i) * vtNa;
    JtK0 = param.GtK(i) * vtK;
    JNkcc10 = param.aNkcc1(i) * ( ( param.a1 - param.a2 * Na0 * K0 * Cl0^2 ) ...
                            / ( param.a3 + param.a4 * Na0 * K0 * Cl0 )^2 );
    vNaK = param.Sb(i) * ( param.r * param.Ke^2 * Na0^3 ...
                           / ( param.Ke^2 + param.alpha1 * Na0^3 ) );                  
    param.aNaK(i) = ( ( JtNa0 + JtK0 ) - JNkcc10 ) / ( 3 * vNaK);
    
    
    % Ca2+ activated K+ channel GK
    PKb = 1 ./ ( 1 + ( param.KCaKC ./ c0_init ).^param.eta2 );
    vK = PKb * ( Vb0 - VK0 ) / param.F;
    param.GK(i) = ( JNkcc10 + 2 * ( JtNa0 + JtK0 ) ) / ( 3 * vK );
    
    % Anion exchanger 4 G4
    JNaK0 = param.aNaK(i) * vNaK;
    JNhe10 = param.G1(i) * vNhe1;
    vAe4 = param.Sb(i) *  ( ( param.Cle / ( param.Cle + param.KCl ) ) * ...
        ( Na0 / ( Na0 + param.KNa ) ) ...
             * ( HCO30 / ( HCO30 + param.KB ) )^2 );
    param.G4(i) = ( JNkcc10 - 3 * JNaK0 + JNhe10 ) / vAe4;
    
    param.GtNa(i) = 10*param.GtNa(i);
    param.GtK(i) = 10*param.GtK(i);
    param.GCl(i) = 0.7*param.GCl(i);
    param.aNaK(i) = param.aNaK(i)*1.1;
end
save('single_cell_conductance_file.mat','param','w_IPR','w_basal')



%% Computing
par = cell(n_cells,1);
dt = zeros(n_cells,1);
for i = 1:n_cells
    dt(i)= delt;
    par{i}= param;
end
tic
% parpool(7)

% Inititalizing time
time = 0*dt;

% Tolerance for the variable time step
tol = 1e-3;
t=0;

% Indices for the stored time steps1
ind =zeros(n_cells,1);
time_ind = cell(n_cells,1);
for i = 1:n_cells
    time_ind{i} =time(i);
end

% Initialization of the cells containing calcium, ip3 and ions
x = cell(n_cells,1);
x_RyR = cell(n_cells,1);
x_h = cell(n_cells,1);
x_ion = cell(n_cells,1);
vdt = zeros(n_cells,1);
x_tilde = cell(n_cells,1);


% To know if still computing
keep_going = true(n_cells,1);

% Rigidity matrix cell
Amat2 = cell(n_cells,1);


% Initial values
for i = 1:n_cells
    x{i} = u{i}(:,1);
    x_RyR{i} = h_RyR{i}(:,1);
    x_tilde{i} = x{i};
    x_h{i} = h_border{i}(:,1);
    x_ion{i} = ion{i}(:,1);
    vdt(i) = volume(i);
end

while ((min(ind)<(numt-1))&& (min(time(:,end))<tend))
    t = t+1;
    time = [time,time(:,end)+dt];

    for i = 1:n_cells
        Amat{i} = mass{i}+dt(i)*stiff{i};
        Amat2{i} = mass{i}+(dt(i)/2)*stiff{i};
    end
    e_mean = 0;
    e_max = 0;

% To run in parallel over each cell
    for i = 1:n_cells
        if time(i,end)<tend
            
% Set VPLC=0 for the first few seconds
        if (time(i,end) < param.initialperiod)
            param.VPLC = 0;
            w_VPLC{i} = param.VPLC*(dist_ba{i}< param.PLCds).*(dist_ap{i}> param.PLCdl);
        else
            param.VPLC = param.VPLCkeep;
            w_VPLC{i} = param.VPLC*(dist_ba{i}< param.PLCds).*(dist_ap{i}> param.PLCdl);
        end

% Ion
        c=x_tilde{i}(1:np(i));      
        cav_tri  = m_triangles{i}(:,1:size(x_h{i},1))*c(1:size(x_h{i},1));
         
%  Here is the first step that couples the calcium model to secretion        
        f_secretion = @(t,x) Secretion_new(t,x,cav_tri,param, surf_tri{i},i,...
            w_IPR{i},w_basal{i});
        [~,ion_fstep] = ode15s(f_secretion,[0 dt(i)/2],x_ion{i}'); 
        
        vdt2 = ion_fstep(end,4);
%             e = max(abs(w-x_ion{i}));
        

        % Volume 
        % reaction for the first propagation

        [reac1,reac1_h,reac1_RyR] = make_load(x_tilde{i},x_h{i}, x_RyR{i}, p{i},...
            par{i},w_IPR{i}, w_VPLC{i}, w_RyR{i}, vol_tets{i},surf_tri{i},...
            m_tets{i},m_triangles{i},volume(i)/vdt2);


        % Main step
        w = Amat{i}\...
            (mass{i}*x{i}+dt(i)*reac1);
%             w_tilde = w;
%             w_tilde(1:2*np(i)) = (volume(i)/vdt(i))*w(1:2*np(i));


        % First half step
        v = Amat2{i}\...
            (mass{i}*x{i}+(dt(i)/2)*reac1);

        v_h = x_h{i} + (dt(i)/2)*reac1_h; 
        v_hRyR = x_RyR{i} + (dt(i)/2)*reac1_RyR;
        v_tilde =v;
        v_tilde(1:2*np(i)) = (volume(i)/vdt2)*v(1:2*np(i));

%  Here is the second step that couples the calcium model to secretion        

        cav_tri  = m_triangles{i}(:,1:size(x_h{i},1))*v_tilde(1:size(x_h{i},1));
        f_secretion = @(t,x) Secretion_new(t,x,cav_tri,param, surf_tri{i},i,...
            w_IPR{i},w_basal{i});
        [~,ion_sstep] = ode15s(f_secretion,[0 dt(i)/2],ion_fstep(end,:)'); 
        vdt(i) = ion_sstep(end,4);
        x_ion{i} = ion_sstep(end,:)';

        % reaction for the second propagation
        [reac,reac_h,reac_RyR] = make_load(v_tilde,v_h, v_hRyR, p{i},...
            par{i},w_IPR{i}, w_VPLC{i}, w_RyR{i},vol_tets{i},surf_tri{i},...
            m_tets{i},m_triangles{i},volume(i)/vdt(i));

        % Second half step
        x{i} = Amat2{i}\...
            (mass{i}*v+(dt(i)/2)*reac);
        x_h{i} = v_h + (dt(i)/2)*reac_h;
        x_RyR{i} = x_RyR{i} + (dt(i)/2)*reac_RyR;
        x_tilde{i} = x{i};
        x_tilde{i}(1:2*np(i)) = (volume(i)/vdt(i))*x{i}(1:2*np(i));

        % errors
        e = max(max(abs(w-x{i})));
        e_mean = e_mean + sum(abs(w-x{i}));
        e_max = max(e,e_max);
        disp(['dt = '...
            num2str(dt(i)),', time t = ',num2str(time(i,end))])
        dt(i) = min(max(min(0.9*sqrt(tol/e)*dt(i),2*dt(i)),0.001),1);
        end
    end

    for i = 1:n_cells
    if (keep_going(i))&& ((time(i,end)-time_ind{i}(end))>=0.1)...
            &&(time(i,end)>(tend-t_keep))
        ind(i) = ind(i) + 1;       
        u{i}(:,ind(i)) = x_tilde{i};
        ion{i}(:,ind(i)) = x_ion{i};
        time_ind{i} = [time_ind{i},time(i,end)];
        %disp(['keep cell ',num2str(i)])
    end
        
        
    end

e_max = max(e,e_max);
    
%% 

for i = 1:n_cells
    if (time(i,end)>tend)
        if keep_going(i)
            keep_going(i) = false;
            c_tot{i} = u{i}(1:np(i),1:size(time_ind{i},2));
            c_tot{i} = interp1(time_ind{i},c_tot{i}',reg_time)';
            ip_tot{i} = u{i}(np(i)+1:2*np(i),1:size(time_ind{i},2));
            ip_tot{i} = interp1(time_ind{i},ip_tot{i}',reg_time)';
            ion_tot{i} = interp1(time_ind{i},...
                ion{i}(:,1:size(time_ind{i},2))',reg_time)';
            disp('interpolating')
            u{i} =0;
            x{i} = 0;
            x_RyR{i} = 0;
            x_h{i} = 0;
            vdt(i) = 0;
        end
    end
end

%disp(['Step n ',num2str(t),', mean error =',num2str(e_mean/sum(np)),...
%        ',max error  = ',num2str(e_max),', elapsed time ', num2str(toc)])
end
toc



%% plotting results

    save('single_cell_output','reg_time','c_tot','ip_tot','ion_tot','-v7.3')


