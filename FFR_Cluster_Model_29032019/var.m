function [int,Nal,Kl,Cll] = var(x,par)
% #########################################################################
% -------------------------------------------------------------------------
% Author: Elias Siguenza
% Location: The University of Auckland, New Zealand
% Date: 22 March 2019
% Version: 2.1
% -------------------------------------------------------------------------
% Purpose:
% This function sorts out the variables to a matrix of seven columns
% (that represent each cell) and 8 rows (representing concentrations and PM
% potentials along with cell volume). 
% The order is: 
% Volume, 
% Na+ 
% K+
% Cl- 
% HCO3- 
% H+
% Va 
% Vb
% -------------------------------------------------------------------------
% #########################################################################
% Initialise and preallocate memory for the loop.
int = zeros(8,7);
Nal = zeros(7,7);
Kl = zeros(7,7);
Cll = zeros(7,7);
n = 1;
l = 57;
%%%%%%%%%%%%%%%%% Begin Loop
for i = 1:7
    int(1:8,i) = x(n:n+7);
    n= n+8;
    for j = 1: par.ind_clust{i}
        ngh = par.neigh_clust{i}(j);
            Nal(i,ngh) = x(l);
            Kl(i,ngh) = x(l+sum(cell2mat(par.ind_clust)));
            Cll(i,ngh) = x(l+(2*sum(cell2mat(par.ind_clust))));
            l = l+1;
    end
end
%%%%%%%%%%%%%%%%%
% Luminal concentration matrix:
% |1,0,0,0,0,0,0|
% |1,1,0,0,0,0,0|
% |1,1,1,0,0,0,0|
% |1,1,1,1,0,0,0|
% |0,1,1,1,0,0,0|
% |0,1,0,1,1,1,0|
% |0,0,0,0,1,1,0|
%%%%%%%%%%%%%%%%%
% NOTE: The luminal concentration matrix is
% essentially a lower triangular matrix whose rows represent cell number
% and its columns represent neighbour. However, for  calculation of
% membrane potentials and tight junctional fluxes some of these must be
% repeated. However they are not variables of the system. To get around
% that, I just add the matrix's transpose to the variable matrix and I
% obtain what I want:

Nal = Nal + (tril(Nal)'-diag(diag(Nal)));
Kl = Kl + (tril(Kl)'-diag(diag(Kl)));
Cll = Cll + (tril(Cll)'-diag(diag(Cll)));
%%%%%%%%%%%%%%%%%

end