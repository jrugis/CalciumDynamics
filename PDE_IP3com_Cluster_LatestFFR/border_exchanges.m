function [exchange_reaction] = border_exchanges(x_p,m_triangles,...
                            com_tri,param,np,nptri,surf_tri)
% This function implements the exchanges of calcium and IP_3 between the
% cells of the cluster

% Input:
% x_p : cell of vector of the concentrations inside the cell, first np(i)
% are calcium concentration inside the cell, from np(i)+1 to 2*np(i) IP_3
% concentration and from 2*np(1)+1 to 3*np(i) calcium concentration in the
% ER.
% m_triangles : cell of matrices that make the transformation from calcium
% concentration at each node to each triangle.
% com_tri: n_cell by n_cell cell that include each of the triangles shared
% between cell i and j. by indices
% triangles : cell of vectors that include the indices of the nodes for
% each triangle
% param: structure including the parameters
% np : vector of number of nodes in each cell
% nptri : vectornumber of points that are on the border

% Output:
% exchange_reaction : cell of vector including the reaction for each cell

n_cells = size(x_p,1);
c = cell(n_cells,1);
ip = cell(n_cells,1);

cav_tri = cell(n_cells,1);
ipav_tri = cell(n_cells,1);

exchange_reaction = cell(n_cells,1);
load_ip = cell(n_cells,1);
for i = 1:n_cells
    
% Redifining the input so it is more conveniejnt
%     c{i}=u(1:np(i));
    ip{i}=x_p{i}(np(i)+1:2*np(i));


% Computing the values on the triangles
%     cav_tri{i}  = m_triangles{i}(:,1:nptri(i))*c{i}(1:nptri(i));
    ipav_tri{i} = m_triangles{i}*ip{i};
    
% Initialise exchanges
exchange_reaction{i} = zeros(2*np(i),1);
load_ip{i} = zeros(np(i),1);
end



for i = 1:n_cells
    for j = i+1:n_cells
        if ~isempty(com_tri{i,j})

            ip_tri_exchanges = param.Fip * (ipav_tri{j}(com_tri{i,j}(:,2))-...
                ipav_tri{i}(com_tri{i,j}(:,1)));
            

            ip_tri_i = zeros(nptri(i),1);
            ip_tri_i(com_tri{i,j}(:,1)) = ip_tri_exchanges;

            ip_tri_j = zeros(nptri(j),1);
            ip_tri_j(com_tri{i,j}(:,2)) = -ip_tri_exchanges;
            

            load_ip{i} = load_ip{i} + (m_triangles{i}')*(surf_tri{i}.*ip_tri_i);
            load_ip{j} = load_ip{j} + (m_triangles{j}')*(surf_tri{j}.*ip_tri_j);
        end
    end
end

for i = 1:n_cells
    exchange_reaction{i}(np(i)+1:2*np(i)) = load_ip{i};
end

