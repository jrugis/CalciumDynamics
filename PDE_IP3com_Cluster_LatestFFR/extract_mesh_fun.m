function [p, tets, volume, tets_volume, dup_nodes,tri, bndry,dist_apical,dist_ap_p,dist_basal] = ...
    extract_mesh_fun(files,shw, n_cells_max, mod_basal)
% The goal of this script is to extract the mesh data form the ASCII file.
% Then to duplicate the nodes at the fronter between two cells and to
% rearrange the tetrahedrons with the new notations. The end of the script
% plot the different cells.

%% Input 


%% Output


%% Reading file

% extracting raw information from each of the cell mesh
n_cells =min(size(files,2),n_cells_max);
np = cell(n_cells,1);
n_tri = cell(n_cells,1);
n_tets = cell(n_cells,1);
p_ind = cell(n_cells,1);
p = cell(n_cells,1);
tri = cell(n_cells,1);
tri_ind = cell(n_cells,1);
tets = cell(n_cells,1);
tets_ind = cell(n_cells,1);
dist_ap_p = cell(n_cells,1);
dist_apical = cell(n_cells,1);
for i =1:n_cells
    disp(['Extracting mesh infos for cell ',num2str(i)])
    np{i} = dlmread(files{i}, ' ', [4 0 4 0]);
    % reading the nodes
    p_ind{i} = dlmread(files{i}, ' ', [5 0 5+(np{i}-1) 0]);
    % reading to which cell the nodes belong to
    p{i} = dlmread(files{i}, ' ', [5 1 5+(np{i}-1) 3]);
    
    n_elem =  dlmread(files{i}, ' ', [5+(np{i}-1)+3 0 ...
                            5+(np{i}-1)+3 0]);
    tri_tet = dlmread(files{i}, ' ', [5+(np{i}-1)+4 0 ...
                            5+(np{i}-1)+4+(n_elem-1) 8]);
    n_tri{i} = find(tri_tet(:,5)~=tri_tet(1,5),1)-1;
    tri{i} = tri_tet(1:n_tri{i},6:8);
    tri_ind{i} = tri_tet(1:n_tri{i},1);
    
    n_tets{i}= size(tri_tet,1)-n_tri{i};
    % reading the nodes composing the tetrahedron
    tets{i} = tri_tet(n_tri{i}+(1:n_tets{i}),6:9);
    % reading to which cell belong the tets 
    tets_ind{i} = tri_tet(n_tri{i}+(1:n_tets{i}),1);
    
    
    % Distance to apical is now computed using the simplified mesh
%     dist_ap_p{i} = dlmread(files{i}, ' ', [5+(np{i}-1)+4+(n_elem-1)+11 1 ...
%                             5+(np{i}-1)+4+(n_elem-1)+11+np{i}-1 1]);
%     dist_apical{i} = sum(dist_ap_p{i}(tets{i}),2)/4;
end

% Exctrating information from simplified lumen
tree_nodes = dlmread('tree.txt', ' ', [1 0 73 2]);
segments = dlmread('tree.txt', ' ', [75 0 146 1]);
                    
clear n_elem
clear tri_tet 



%% Computing the distance to the apical region
disp('Computing apical distances')
for i = 1:n_cells
    disp(['For cell ',num2str(i)])
    dist_ap_p{i} = 10 * ones(np{i},1);
    for seg = 1:size(segments,1)
            vector_segment = tree_nodes(segments(seg,2),:) - ...
                tree_nodes(segments(seg,1),:);
        for node =1:np{i}
            vector_1 = p{i}(node,:)-tree_nodes(segments(seg,1),:);
            vector_2 = p{i}(node,:)-tree_nodes(segments(seg,2),:);
            % verify using scalar product that the point is aligned with
            % the segment
            scal_1 = sum(vector_segment.*vector_1);
            scal_2 = -sum(vector_segment.*vector_2);
            
            if ((scal_1>0)& (scal_2>0))           
                dist = norm(cross(vector_segment,vector_1))...
                    /norm(vector_segment);
                dist_ap_p{i}(node) = min(dist_ap_p{i}(node),dist);
            elseif ((scal_1>0)& (scal_2<0))
                dist = norm(vector_2);
                dist_ap_p{i}(node) = min(dist_ap_p{i}(node),dist);
            elseif ((scal_1<0)& (scal_2>0)) 
                dist = norm(vector_1);
                dist_ap_p{i}(node) = min(dist_ap_p{i}(node),dist);
            else
                disp('Something is wrong in the distance to apical region')
            end
        end
    end
    dist_apical{i} = sum(dist_ap_p{i}(tets{i}),2)/4;
end

%% Duplicating nodes at borders


% This matrix will contain the number of common nodes between cell i and j.
% (i,i) cells are irrelevant
cmn_bndry = zeros(n_cells);
cmn_p = cell(n_cells);

for i = 1:n_cells
    for j = (i+1):n_cells
        disp(['Computing the common nodes between cell ',num2str(i),' and ',num2str(j)])
        dist_nodes = abs(repmat(p{i}(:,1),1,np{j})-...
            repmat((p{j}(:,1))',np{i},1));
        dist_nodes = dist_nodes + abs(repmat(p{i}(:,2),1,np{j})-...
            repmat((p{j}(:,2))',np{i},1));
        dist_nodes = dist_nodes + abs(repmat(p{i}(:,3),1,np{j})-...
            repmat((p{j}(:,3))',np{i},1));
        [cmn_p{i,j},cmn_p{j,i}] = find((dist_nodes<1e-6));
        cmn_bndry(i,j) = size(cmn_p{i,j},1); 
        
    end 
end
cmn_bndry = cmn_bndry+ cmn_bndry';
clear dist_nodes_1 dist_nodes_2 dist_nodes_3 dist_nodes

%% Verification of the border

for i = 1:n_cells
    disp(['Verification of the common nodes for cell ', num2str(i)])
    cmn_nodes = [];
    for j = 1:n_cells
        cmn_nodes = [cmn_nodes;cmn_p{i,j}];
    end
    nb_cmn = sum(ismember(unique(cmn_nodes),tri{i}));
    if (nb_cmn== max(size(unique(cmn_nodes))))
        disp('Ok')
    end
end
%% Detetcting basal region and creating the Matrix of flux boundaries
dup_nodes = cell(n_cells);
for i = 1:n_cells
    for j = (i+1):n_cells
        dup_nodes{i,j} = zeros(np{i}+np{j});
    end
end


bndry = cell(n_cells,1);
% Detects all the triangles at the bndry
for i = 1:n_cells
    disp(['Detecting the external boundaries and the',...
        ' transmission matrix for cell ', num2str(i)])
    bndry{i} = [];
    for k = 1:n_tri{i}
        bol_bndry = true;
        t = tri{i}(k,:);
        for j = 1:n_cells
            if cmn_bndry(i,j)>0
                if (sum(ismember(t,cmn_p{i,j})*1)==3)
                    bol_bndry = false;
                    vert = [p{i}(t(1),:); p{i}(t(2),:); p{i}(t(3),:)];  % Coords of the vertices
                    side1=vert(1,:)-vert(2,:); 
                    side2=vert(1,:)-vert(3,:);
                    triangle_area = 0.5*norm(cross(side1,side2));
                    v = [];
                    for l = 1:3
                        v = [v;find(cmn_p{i,j}==t(l))];
                    end
                    ind_i = min(i,j);
                    ind_j = max(i,j);
                    v1 = cmn_p{ind_i,ind_j}(v);
                    v2 = np{ind_i}+cmn_p{ind_j,ind_i}(v);
                    dup_nodes{ind_i,ind_j}([v1,v2],[v1,v2]) ...
                        = dup_nodes{ind_i,ind_j}([v1,v2],[v1,v2]) + ...
                        (1/3)*triangle_area*[eye(3),-eye(3);...
                                            -eye(3),eye(3)];
                end
            end
        end
        if bol_bndry
            bndry{i} = [bndry{i};t];
        end
    end
end

for i = 1:n_cells
    for j = (i+1):n_cells
        dup_nodes{i,j} = sparse(dup_nodes{i,j}/2);
        if cmn_bndry(i,j)>0
%             figure,
%             spy(dup_nodes{i,j})
%             title([num2str(i),num2str(j)])
            sum(dup_nodes{i,j}(:))
        end
    end
end


%% Compute the volume

disp('Compute the volume')
volume = zeros(n_cells_max,1);
tets_volume = cell(n_cells_max,1);
for i = 1:n_cells_max
    for t = 1:n_tets{i}
        u = p{i}(tets{i}(t,1),:)-p{i}(tets{i}(t,4),:);
        v = p{i}(tets{i}(t,2),:)-p{i}(tets{i}(t,4),:);
        w = p{i}(tets{i}(t,3),:)-p{i}(tets{i}(t,4),:);
        tets_volume{i}(t) = (1/6)*abs(det([u;v;w]));
        volume(i)= volume(i)+tets_volume{i}(t);
    end
end

for i = 1:n_cells_max
    disp(['number of nodes in cell ', num2str(i), ' : ', num2str(size(p{i},1)),' with volume ',num2str(volume(i)) ])    
end

for i = 1:(n_cells_max-1)
    for j = i+1:n_cells
        disp(['number of common nodes between cell ', num2str(i), ' and ', num2str(j), ' : ', num2str(cmn_bndry(i,j))])
    end
end


%% Distance to basal region
dist_basal = cell(n_cells,1);
for i =1:n_cells
    disp(['Computing distance to basal membrane for cell ',num2str(i)])
    if mod_basal
        basal_membrane = unique(bndry{i}(:));
    else
        basal_membrane = unique(tri{i}(:));
    end
    np_basal = max(size(basal_membrane));
    dist_membrane = (repmat(p{i}(:,1),1,np_basal)-...
            repmat((p{i}(basal_membrane,1))',np{i},1)).^2;
    dist_membrane = dist_membrane + (repmat(p{i}(:,2),1,np_basal)-...
            repmat((p{i}(basal_membrane,2))',np{i},1)).^2;
    dist_membrane = dist_membrane + (repmat(p{i}(:,3),1,np_basal)-...
            repmat((p{i}(basal_membrane,3))',np{i},1)).^2;
    dist_basal_p = sqrt(min(dist_membrane,[],2));
    dist_basal{i} = mean(dist_basal_p(tets{i}),2);
end

clear dist_membrane dist_basal_p
%% Plot 

if shw


    for i = 1:n_cells_max
        
        figure,
        trisurf(tri{i},p{i}(:,1),p{i}(:,2),p{i}(:,3))
        title(['Cell number : ', num2str(i)])
    end

end





end

