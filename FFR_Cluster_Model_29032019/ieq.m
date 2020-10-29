function eq = ieq(JCL,JtNa,JtK,Jb,int,c_no)
% #########################################################################
% -------------------------------------------------------------------------
% Author: Elias Siguenza
% Location: The University of Auckland, New Zealand
% Date: 22 March 2019
% Version: 2.1
% -------------------------------------------------------------------------
% Purpose:
% This function calculates the differential equations of the intracellular
% ionic concentrations and the membrane potentials of any given cell
% using the int matrix of seven columns
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
% dw/dt
eq(1,1) = Jb(c_no,9);
% d[Na+]i/dt
eq(2,1) = (Jb(c_no,3)-3*Jb(c_no,2)+Jb(c_no,5)-Jb(c_no,4)-Jb(c_no,9)*int(2,c_no))/int(1,c_no);
% d[K+]i/dt
eq(3,1) = (Jb(c_no,3)+2*Jb(c_no,2)-Jb(c_no,7)-Jb(c_no,9)*int(3,c_no))/int(1,c_no);
% d[Cl-]i/dt
eq(4,1) = (2*Jb(c_no,3)+Jb(c_no,4)+JCL(c_no,1)-Jb(c_no,9)*int(4,c_no))/int(1,c_no);
% d[HCO3-]i/dt
eq(5,1) = (Jb(c_no,6)-2*Jb(c_no,4)-Jb(c_no,9)*int(5,c_no))/int(1,c_no);
% d[H+]i/dt
eq(6,1) = (Jb(c_no,6)-Jb(c_no,5)-Jb(c_no,9)*int(6,c_no))/int(1,c_no);
% dVa/dt
eq(7,1) = - JCL(c_no,1) - (sum(JtNa(c_no,:)) + sum(JtK(c_no,:)));
% dVb/dt
eq(8,1) = - Jb(c_no,2) - Jb(c_no,7) + (sum(JtNa(c_no,:)) + sum(JtK(c_no,:)));
end