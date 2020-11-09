function [mass, stiff] = make_matrices(p, tets, Dc, Dp)

n_cell = size(p,1);
n_tets = zeros(n_cell,1);
np = zeros(n_cell,1);
for i= 1:n_cell
    np(i) = size(p{i},1);
    n_tets(i) = size(tets{i},1);
end
stiff = cell(n_cell,1);
mass = cell(n_cell,1);

for i = 1:n_cell
    disp(['Making stiffness and mass matrix for cell ', num2str(i)])
    
    stiffc=spalloc(np(i),np(i),8*np(i));
    stiffp=spalloc(np(i),np(i),8*np(i));
    small_mass=spalloc(np(i),np(i),8*np(i));
    refmass=zeros(4);

    % Make the reference mass matrix by integrating over the reference
    % tetrahedron. Just for now, do a totally kludgy integral as a one-point
    % formula. Nasty, and should be improved.

     for j=1:4
        for k=1:4
            ph1=str2func(strcat('phi_',num2str(j)));
            ph2=str2func(strcat('phi_',num2str(k)));
            refmass(j,k)=(1/6)*ph1(0.25,0.25,0.25)*ph2(0.25,0.25,0.25);
        end
     end

     %dum=sum(refmass,2);
     %refmass=diag(dum,0);

    % Now go over every triangle, making the stiffness, load and mass matrices.
    for k=1:n_tets(i)
        vi=tets{i}(k,1:4);    % These are the numbers of the vertices of the tetrahedral element.
        vert = [p{i}(vi(1),:); p{i}(vi(2),:); p{i}(vi(3),:); p{i}(vi(4),:)];  % Coords of the vertices
        J=[vert(2,:)'-vert(1,:)' vert(3,:)'-vert(1,:)'  vert(4,:)'-vert(1,:)'];

        M = [ones(4,1)  vert];
        C = inv(M);
        G=C(2:4,:)'*C(2:4,:);  % Gradients of the basis functions, stored in a matrix

        detJ=abs(det(J));
        A=detJ/6;  % Area of Tet k
        Ic=A*Dc;  % This assumes that the diffusion coefficient is just constant
        Ip=A*Dp;

        % Now put the components of G*I into the correct bits of the stiffness
        % matrix. I can't be bothered writing loops to do this 
        stiffc(vi,vi)=stiffc(vi,vi)+G*Ic;
        stiffp(vi,vi)=stiffp(vi,vi)+G*Ip;


    % Now construct the global mass matrix


        small_mass(vi,vi)=small_mass(vi,vi)+refmass*detJ;

    end  %  of the k loop over the tetrahedrons
    

    % For every variable that doesn't diffuse, set the corresponding elements
    % of the mass matrix just to be the identity matrix
    mass{i} = blkdiag(small_mass,small_mass);
    stiff{i} = blkdiag(stiffc,stiffp);
end


end


%-----------------------------------------

function out=phi_1(s,t,u)
out=1.0-s-t-u;
end
%-----------------------------------------

function out=phi_2(s,t,u)
out=s;
end
%-----------------------------------------

function out=phi_3(s,t,u)
out=t;
end
%-----------------------------------------

function out=phi_4(s,t,u)
out=u;
end
