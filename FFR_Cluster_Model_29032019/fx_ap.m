function [JtNa,JtK,JCl,Qa, Qtot] = fx_ap(int,Nal,Kl,Cll,Jb,Ca,par,c_no,ngh) 
% #########################################################################
% -------------------------------------------------------------------------
% Author: Elias Siguenza
% Location: The University of Auckland, New Zealand
% Date: 22 March 2019
% Version: 2.1
% -------------------------------------------------------------------------
% Purpose:
% This function takes as an input the intracellular and luminal variables 
% of the cell model and the Ca in order to calculate the apical 
% mebrane ionic fluxes, and the flow rate out and into the lumen, 
 % of any particular cell. 
% -------------------------------------------------------------------------
% #########################################################################
% Tight Junctional Membrane Potential
Vt = int(7,c_no)-int(8,c_no);
% -------------------------------------------------------------------------
% Ca2+ Activated Apical Cl- Channels
PCl=1./(1+(par.KCaCC./Ca{1}{c_no,ngh}).^par.eta1);
PrCl = sum(PCl.*par.com_tri_ap{c_no,ngh}(:,3))/par.Sa{c_no};
VCl = par.Sa_p{c_no,ngh} * par.RTF * log(Cll(c_no,ngh)/int(4,c_no));
JCl= par.GCl * PrCl * (int(7,c_no) + (VCl/par.Sa_p{c_no,ngh}))/par.F;
% -------------------------------------------------------------------------
% Tight Junctional Fluxes
VtNa = par.Sa_p{c_no,ngh} * par.RTF * log(Nal(c_no,ngh)/par.Nae);
JtNa = par.Sa_p{c_no,ngh}*par.GtNa*par.St*(Vt-(VtNa/par.Sa_p{c_no,ngh}))/par.F;
% -------------------------------------------------------------------------
VtK = par.Sa_p{c_no,ngh} * par.RTF * log(Kl(c_no,ngh)/par.Ke);
JtK = par.Sa_p{c_no,ngh}*par.GtK*par.St*(Vt-(VtK/par.Sa_p{c_no,ngh}))/par.F;
% -------------------------------------------------------------------------
% Luminal Osmolarity (Using electroneutrality principle)
Il = 2*Cll(c_no,ngh) + par.Ul;
% -------------------------------------------------------------------------
% Flow Rates
% Apical Flow Rate
Qa = par.Sa_p{c_no,ngh} * par.La * (Il - Jb);
% Tight Junctional Flow Rate
Qt = par.Sa_p{c_no,ngh} * par.Lt * (Il - par.Ie);
% Total Fluid Flow Rate (into/out of lumen)
Qtot =(Qa + Qt);
% #########################################################################
end