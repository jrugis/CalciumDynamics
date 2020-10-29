function [Qb, JNaK, JNkcc1, JAe4, JNhe1, JBB, JK, Ii] = fx_ba(int,Ca,par,cell_no)
% #########################################################################
% -------------------------------------------------------------------------
% Author: Elias Siguenza
% Location: The University of Auckland, New Zealand
% Date: 22 March 2019
% Version: 2.1
% -------------------------------------------------------------------------
% Purpose:
% This function takes as an input the intracellular variables of the cell
% model and the Ca in order to calculate the mebrane ionic fluxes, 
% and the flow rate into the cell, at the basolateral side of any 
% particular cell. 
% -------------------------------------------------------------------------
% #########################################################################
% Basolateral Flow Rate
    Qb = par.Lb * ( 2 * ( int(2,cell_no) + int(3,cell_no) + int(6,cell_no) ) + par.CO20 - par.Ie);
% K+ Nernst Potential    
    VK = par.RTF * log(par.Ke/int(3,cell_no));
% Ca2+ Activated K+ Channels open probability    
    PK = sum((1./(1+(par.KCaKC./(Ca{2}{cell_no})).^par.eta2)).*par.Sb_k{cell_no})./par.Sb{cell_no};
% Ca2+ Activated K+ Channels Total Flux    
    JK = par.GK * PK * ( int(8,cell_no) - VK ) / par.F;
% 3Na+/2K+ ATPases     
    JNaK  = par.Sb{cell_no}*par.aNaK*(par.r*par.Ke^2*int(2,cell_no)^3/(par.Ke^2+par.alpha1*int(2,cell_no)^3));
% Na+ K+ 2Cl- Cotransporters
    JNkcc1= par.aNkcc1*par.Sb{cell_no}*(par.a1-par.a2*int(2,cell_no)*int(3,cell_no)*int(4,cell_no)^2)/(par.a3+par.a4*int(2,cell_no)*int(3,cell_no)*int(4,cell_no)^2);
% Na+-2HCO3-/Cl- Anion Exchanger 4    
    JAe4  = par.Sb{cell_no}*par.G4 *((par.Cle/(par.Cle+par.KCl))*(int(2,cell_no)/(int(2,cell_no)+par.KNa))*(int(5,cell_no)/(int(5,cell_no)+par.KB))^2);
% Na+/H+ Antiporter 1
    JNhe1 = par.Sb{cell_no}*par.G1*((par.Nae/(par.Nae+par.KNa))*(int(6,cell_no)/(par.KH+int(6,cell_no)))-(int(2,cell_no)/(int(2,cell_no)+par.KNa))*(par.He/(par.KH+par.He))); 
% CAIV Intracellular Carbonic Anhydrases    
    JBB   = int(1,cell_no)*par.GB*(par.kp*par.CO20-par.kn*int(5,cell_no)*int(6,cell_no));
% Intracellular Osmolarity (Using electroneutrality principle)    
    Ii = 2*(int(2,cell_no) + int(3,cell_no) + int(6,cell_no)) + par.CO20;
% #########################################################################
end