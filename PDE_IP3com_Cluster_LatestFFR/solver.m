function solver(file_name,folder)

close all
clc

% Loading the parameter file

folder

if folder
	folder_file_name = [num2str(folder),'/',file_name];
	load(folder_file_name)
else
	load(file_name)
end

%% Extract mesh


% Use the function extract_mesh_fun to transform the original mesh and get
% the interesting datas
str = ['mod_basal',num2str(mod_basal),'data_smoothed_mesh.mat'];
if exist(str, 'file') ~= 2
    [p, tets, volume, tets_volume, membranes,triangles, bndry,dist_ap,dist_ap_p,dist_ba] = ...
        extract_mesh_fun(file_mesh, show, n_cells,mod_basal);
    save(str,'p', 'tets', 'volume', 'tets_volume',...
        'membranes','triangles', 'bndry','dist_ap','dist_ap_p','dist_ba')
else
    load(str)
end

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





%% Loading matrices
str = 'matrixFc.mat';
if exist(str, 'file') ~= 2
    [mass, stiff] = make_matrices(p, tets, Dc, Dp);

    disp('Computing the first step matrix')
    Amat = cell(n_cells,1);
    for i = 1:n_cells
        Amat{i} = mass{i};
    end
    figure,
    spy(Amat{1})

    sMem = cell(n_cells,n_cells,2);
    sMass = cell(n_cells);
    for i =1:n_cells
        for j = (i+1):n_cells
            sMass{i,j} = blkdiag(mass{i}(1:np(i),1:np(i)),mass{j}(1:np(j),1:np(j)));
        end
    end


    save(str,'Amat','sMass','sMem','mass','stiff')
else
    load(str)
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
ion = zeros(length(IC),numt);
ion(:,1) = IC;

% Initialisation of the h border variable
h_border = cell(n_cells,1);
for i = 1:n_cells
    h_border{i} = zeros(size(unique(triangles{i}(:)),1),numt);
    h_border{i}(:,1) =  h0(i);
end

% To store the results
c_tot = cell(n_cells,1); 
ip_tot = cell(n_cells,1);
ion_tot = zeros(length(IC),numt+1);
reg_time = (tend-t_keep):0.1:tend;


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
% x_ion = zeros(length(IC),1);
vdt = zeros(n_cells,1);
x_tilde = cell(n_cells,1);


% To know if still computing
keep_going = true(n_cells,1);



x_ion = ion(:,1);
% Initial values
for i = 1:n_cells
    x{i} = u{i}(:,1);
    x_RyR{i} = h_RyR{i}(:,1);
    x_tilde{i} = x{i};
    x_h{i} = h_border{i}(:,1);
    vdt(i) = volume(i);
end

ca_b = cell(7,1);
ca_a = cell(7,7);
%%
for i = 1:n_cells
        Amat{i} = mass{i}+dt(i)*stiff{i};
        nptri(i) = size(m_triangles{i},1);
end


    
while ((min(ind)<(numt-1))&& (min(time(:,end))<tend))
    t = t+1;
    time = [time,time(:,end)+dt];

        for j = 1:n_cells
            c=x_tilde{j}(1:np(j));
            cav_tri  = m_triangles{j}(:,1:size(x_h{j},1))*c(1:size(x_h{j},1));
            Ca(j) = mean(cav_tri);
            ca_b{j} = ones(length(pa.sb_tri{j}),1)*Ca(j);
            for ii = 1:n_cells
                if ~isempty(pa.sa_tri{j,ii})
                    ca_a{j,ii}=ones(length(pa.sa_tri{j,ii}),1)*Ca(j);
                else
                    ca_a{j,ii}=[];
                end
            end
        end
        ca_both = {ca_a, ca_b};
        
 %      f_secretion = @(t,x) Cluster(t,x,ca_a,ca_b,clust_cell,pa);
       f_secretion = @(t,x) cluster(t, x, ca_both, pa);
       [~,ion_fstep] = ode15s(f_secretion,[0 dt(i)],x_ion'); 
       x_ion = ion_fstep(end,:);
        
%         vdt2(7)    = ion_fstep(end,1);
%         vdt2(6)    = ion_fstep(end,7);
%         vdt2(5)    = ion_fstep(end,13);
%         vdt2(4)    = ion_fstep(end,19);
%         vdt2(3)    = ion_fstep(end,25);
%         vdt2(2)    = ion_fstep(end,31);
%         vdt2(1)    = ion_fstep(end,37);
        vdt2(1) = ion_fstep(end, 1);
        vdt2(2) = ion_fstep(end, 9);
        vdt2(3) = ion_fstep(end, 17);
        vdt2(4) = ion_fstep(end, 25);
        vdt2(5) = ion_fstep(end, 33);
        vdt2(6) = ion_fstep(end, 41);
        vdt2(7) = ion_fstep(end, 49);
        
        load_ip = border_exchanges(x_tilde,m_triangles,...
                            com_tri,param,np,nptri,surf_tri) ;  
    for i = 1:n_cells

        [reac1,reac1_h,reac1_RyR] = make_load(x_tilde{i},x_h{i}, x_RyR{i}, p{i},...
            par{i},w_IPR{i}, w_VPLC{i}, w_RyR{i}, vol_tets{i},surf_tri{i},...
            m_tets{i},m_triangles{i},volume(i)/vdt2(i));

        % Main step
        x{i} = Amat{i}\...
            (mass{i}*x{i}+dt(i)*(reac1+load_ip{i}));


        x_h{i} = x_h{i} + (dt(i))*reac1_h; 
        x_RyR{i} = x_RyR{i} + (dt(i)/2)*reac1_RyR;
        x_tilde{i} = x{i};
        x_tilde{i} = (volume(i)/vdt2(i))*x{i};
        

    end   
        disp(['Dt = ', num2str(dt(i)),...
            ', Simulation Time = ',num2str(dt(i)*t)])

    for i = 1:n_cells
    if (keep_going(i))&& ((time(i,end)-time_ind{i}(end))>=0.1)...
            &&(time(i,end)>(tend-t_keep))
        ind(i) = ind(i) + 1;       
        u{i}(:,ind(i)) = x_tilde{i};
        ion(:,ind(i)) = x_ion;
        time_ind{i} = [time_ind{i},time(i,end)];
        disp(['keep cell ',num2str(i)])
    end
  
    end


for i = 1:n_cells
    if (time(i,end)>tend)
        if keep_going(i)
            keep_going(i) = false;
            c_tot{i} = u{i}(1:np(i),1:size(time_ind{i},2));
            c_tot{i} = interp1(time_ind{i},c_tot{i}',reg_time)';
            ip_tot{i} = u{i}(np(i)+1:2*np(i),1:size(time_ind{i},2));
            ip_tot{i} = interp1(time_ind{i},ip_tot{i}',reg_time)';
            disp('interpolating')
            u{i} =0;
            x{i} = 0;
            x_RyR{i} = 0;
            x_h{i} = 0;
            vdt(i) = 0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (time(i,end)>tend)
    ion_tot = interp1(time_ind{1},ion(:,1:size(time_ind{1},2))',reg_time)';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['Step n ',num2str(t),', elapsed time ', num2str(toc)])
end
toc


%% plotting results
if folder
	sol_folder_file_name = [num2str(folder),'/sol',file_name];
else
	sol_folder_file_name = ['sol',file_name]; 
end
save(sol_folder_file_name,'reg_time','c_tot','ip_tot','ion_tot','-v7.3')

end

