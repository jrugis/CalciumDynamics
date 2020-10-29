function dx = cluster(~,x,Ca,par)
% #########################################################################
% -------------------------------------------------------------------------
% Author: Elias Siguenza
% Location: The University of Auckland, New Zealand
% Date: 22 March 2019
% Version: 4.1
% -------------------------------------------------------------------------
% Purpose:
% This function solves the differential equation system involved in the
% simulation of flow rate and ionic homeostasis in seven parotid acinar 
% cells. The cells are connected according the connectivity matrix 
% determined by the mesh created by John Rugis. 
% This function requires an initial condition and a Calcium input. 
% Such Ca input is given by Nathan Pages' spatial [Ca2+]_i model.
% -------------------------------------------------------------------------
% #########################################################################
%%%% Variables
% Sort out the variables (Arguably this could change, but I like order).
% The order is: 
% 1. Volume 
% 2. Na+ 
% 3. K+ 
% 4. Cl- 
% 5. HCO3- 
% 6. H+ 
% 7. Va 
% 8. Vb
    [int,Nal,Kl,Cll] = var(x,par);
% -------------------------------------------------------------------------
%%%% Memory Preallocation  
    Jb = zeros(7,9);JCl = zeros(7,7);JCL= zeros(7,1);JtNa = JCl;JtK = JCl;
    Qa=JCl; Qtot=JCl;a = 1;dx=zeros(length(x),1);
% -------------------------------------------------------------------------
%%%% Begin calculation of submodels.
    for c_no = 1:7
%%%% Basolateral fluxes and Membrane Potentials
% Sort out the basolateral fluxes:  
% The order is: Qb, JNaK, JNkcc1, JAe4, JNhe1, JBB, JK, Ii, and Jwater.
        [Jb(c_no,1), Jb(c_no,2),...
            Jb(c_no,3), Jb(c_no,4),...
                       Jb(c_no,5), Jb(c_no,6),...
                                    Jb(c_no,7),Jb(c_no,8)] ...
                                                  = fx_ba(int,Ca,par,c_no); 
% -------------------------------------------------------------------------
%%%% Apical fluxes
% Sort out apical and tight junctiuonal fluxes:
        for j = 1:par.ind{c_no}
            ngh = par.neigh{c_no}(j);
            [JtNa(c_no,ngh), JtK(c_no,ngh), JCl(c_no,ngh),...
                                        Qa(c_no,ngh), Qtot(c_no,ngh)]...
                        = fx_ap(int,Nal,Kl,Cll,Jb(c_no,8),Ca,par,c_no,ngh);
        end 
% -------------------------------------------------------------------------        
%%%% Sort out the fluxes        
% Sum all apical chloride fluxes:
        JCL(c_no,1) = sum(JCl(c_no,:));
        Jb(c_no,9) = Jb(c_no,1) - sum(Qa(c_no,:));
% -------------------------------------------------------------------------
%%%% Set the intracellular concentration differential equations
        dx(a:a+7,1) = ieq(JCL,JtNa,JtK,Jb,int,c_no);
        a=a+8;
    end
% -------------------------------------------------------------------------
%%%% Set the luminal concentration equations:
% Use adjacency matrix for connectivity of the lumen.
 [JtNad,JtKd,JCld,Qtotd,Nald,Kld,Clld,QwNa,QwK,QwCl] = ...
                             lum_adj(Nal,Kl,Cll,JtNa,JtK,JCl,Qtot,par);                          
% Initialise loop indices:
    ina = length(int)*7;
    ik = ina+19;
    icl = ik+19;
    
    for i = 1:19
        dx(ina+i,1) = JtNad(i,i) + sum(QwNa(:,i)) ...
                                 -  sum(Qtotd(:,i)) * Nald(i,i); % Sodium
        dx(ik+i,1)  = JtKd(i,i) + sum(QwK(:,i)) ...
                                -  sum(Qtotd(:,i)) * Kld(i,i);  % Potassium
        dx(icl+i,1) = - JCld(i,i) + sum(QwCl(:,i)) ...
                                  -  sum(Qtotd(:,i)) * Clld(i,i);% Chloride
    end
% -------------------------------------------------------------------------    
%%%%
% End of calculations.
end