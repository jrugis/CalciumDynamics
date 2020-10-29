classdef Cell
%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% #### Version 2 ####
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
% Elias Vera-Siguenza
% The University of Auckland
% New Zealand
% 18 November 2018
% Multiscale Modelling the Primary Fluid Secretion
% :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%
    properties
        c_no % Cell Number
        ngh  % Neighbours
        ngh_clust % Neighbours (effective)
        pa   % Parameters
        var  % Variables
        fx   % Fluxes
    end
%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%
    methods % Constructor
        %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
        function c = Cell(cell_number,parameters)
            c.c_no = cell_number;
            c.pa = parameters;
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            c.ngh = c.pa.neigh{cell_number};
            c.ngh_clust = c.pa.neigh_clust{cell_number};
            
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
        end
    end
%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%
% Define the reactions that happen inside the cell:
    methods(Static)
        function [var, fx] = fluxes(c_no,ngh,pa,var_vec,ca_a,ca_b,cluster)
            
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % Sort out the variables according to the connectivity matrix
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % Do you want to compute a cluster, or individual cells?
            % 1 = Cluster, 0 = Single.
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            
            fx.index = pa.ind{c_no};
            
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % This is a bit messy so var_vec must change, 
            % it needs to be automatic.
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            
            var.w = var_vec(1);
            var.Na = var_vec(2);
            var.K = var_vec(3);
            var.Cl = var_vec(4);
            var.HCO3 = var_vec(5);
            var.H = var_vec(6);
            for i = 1:fx.index
                var.Nal(i)=var_vec(6+i);
                var.Kl(i)=var_vec(6+fx.index+i);
                var.Cll(i) = var.Nal(i)+var.Kl(i);
            end
            
            fx.ind_lum = length(var.Nal)*2 + 6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transporters/Pumps/Exchanger in the Basolateral PM:
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            %  Na/K ATPase pump (NaK-ATPase)
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            
            fx.JNaK = pa.Sb{c_no}*pa.aNaK*...
                    (pa.r*pa.Ke^2*var.Na^3/(pa.Ke^2+pa.alpha1*var.Na^3));
            
            % ATPases Current
            iNaK= pa.F * fx.JNaK;
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % Basolateral cotransporter (Nkcc1)
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%

            fx.JNkcc1=pa.aNkcc1 * pa.Sb{c_no} * ...
                  ( pa.a1 - pa.a2 * var.Na * var.K * var.Cl^2 ) ...
                / ( pa.a3 + pa.a4 * var.Na * var.K * var.Cl^2 );
            
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % Anion Exchanger 4 (Ae4)
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            
            fx.JAe4 = pa.Sb{c_no} * pa.G4 ...
                * ( (pa.Cle / (pa.Cle + pa.KCl)) ...
                * (var.Na / (var.Na + pa.KNa)) ...
                * (var.HCO3 / (var.HCO3 + pa.KB))^2);
            
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % Basolateral Na/H+ Antiporter (Nhe1)
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            
            fx.JNhe1 = pa.Sb{c_no} * pa.G1 ...
                    * ( ( pa.Nae / ( pa.Nae + pa.KNa ) ) ...
                            * ( var.H / ( pa.KH + var.H ) ) ...
                                    - ( var.Na / ( var.Na + pa.KNa ) ) ...
                                            *(pa.Hye / ( pa.KH + pa.Hye)));
            
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % Intracellular Carbonic Amylases (This an intracellular
            % reaction)
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            
            fx.JBB = var.w * pa.GB * (pa.kp * pa.CO20 ...
                                               - pa.kn * var.HCO3 * var.H);
            
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % K+ Nernst Potential
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
                
            fx.VK = pa.RTF * log( pa.Ke / var.K );
            
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % Basolateral K+ Channels Open Probability (Mesh-Dependent)
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            
            fx.PrKb=sum((1./(1+(pa.KCaKC./ca_b{c_no}).^pa.eta2)).*pa.Sb_k{c_no})/pa.Sb{c_no};
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Water, Tight Junctional, and Apical PM Fluxes:
%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%            
            % Osmolarities (the interstitial osmolarity is a constant)
            % These are based on electroneutrality.
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % Define the Intracellular Osmolarity Equation
            Ii = 2*(var.Na + var.K + var.H) + pa.CO20;
            % Pre-allocate memory to compute Luminal Osmolarities
            Il = zeros(1,length(ngh));            
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % This computes each flux/osmolarity according to neighbouring
            % cells (given by the connectivity matrix)  
            for j = 1:length(ngh)
                Il(j) = 2*var.Cll(j)+pa.Ul;
                Sa_p = pa.Sa_p{c_no,ngh(j)};
                % Water Fluxes
                fx.Qa(1,j) = Sa_p * pa.La * (Il(j) - Ii);
                %fx.Qt(1,j) = Sa_p * pa.Lt * ( Il(j) - pa.Ie );
                if cluster
                    if c_no >= ngh(j)
                        fx.Qt(1,j) = Sa_p * pa.Lt * ( Il(j) - pa.Ie );
                    else
                        fx.Qt(1,j) = 0; % This is because I don't 
                                        % want to double on the amount of 
                                        % water in the lumen via paracell.
                    end
                else
                    fx.Qt(1,j) = Sa_p * pa.Lt * ( Il(j) - pa.Ie );
                end
                
                fx.Qtot(1,j) = fx.Qa(1,j) + fx.Qt(1,j);
                % Tight Junctional Nernst Potentials
                fx.VtNa(1,j)= Sa_p*pa.RTF*log(var.Nal(j)/pa.Nae);
                fx.VtK(1,j) = Sa_p*pa.RTF*log(var.Kl(j)/pa.Ke );
                fx.VCl(1,j) = Sa_p*pa.RTF* log(var.Cll(j)/var.Cl);
                PCl=1./(1+(pa.KCaCC./ca_a{c_no,ngh(j)}).^pa.eta1);
                fx.PrCl(1,j) = sum( PCl .* ...
                    pa.com_tri_ap{c_no,ngh(j)}(:,3) )/pa.Sa{c_no};
            end

     
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % Basolateral Water Flux.
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            
            fx.Qb = pa.Lb * (Ii - pa.Ie);
            
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % Plasma Membrane Potentials (PM)
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            
            GtX = pa.GtNa * pa.St;
            GtY = pa.GtK * pa.St;
            GY  = pa.GK  * fx.PrKb;
            ND  = pa.St*(pa.GtNa+pa.GtK);
            vtx = sum(fx.VtNa);
            vty = sum(fx.VtK);
            Gz  = pa.GCl * sum(fx.PrCl);
            Gzv = pa.GCl * dot(sum(fx.PrCl),sum(fx.VCl));
            
            fx.Va=-(ND*iNaK-GY*(GtX*vtx+GtY*vty+ND*fx.VK)...
                +Gzv*(GY+ND))/(ND*(GY+Gz)+GY*Gz);
            
            fx.Vb=-((ND+Gz)*(iNaK-GY*fx.VK)+Gz*(GtX*vtx+GtY*vty)...
                +ND*Gzv)/(ND*GY+Gz*(GY+ND));
            Vt = fx.Va - fx.Vb;
            
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % Basolateral Ca2+ activated K Channels
            fx.JK = pa.GK*fx.PrKb*(fx.Vb-fx.VK)/pa.F;
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            
            for i = 1:length(ngh)
                neighbour = ngh(i); Sa_p=pa.Sa_p{c_no,neighbour};
                fx.JCl(1,i) = pa.GCl * fx.PrCl(1,i) * (fx.Va + (fx.VCl(1,i)/Sa_p))/pa.F;
                fx.JtNa(1,i) = Sa_p * pa.GtNa * pa.St * ( Vt - (fx.VtNa(1,i)/Sa_p )) / pa.F;
                fx.JtK(1,i) = Sa_p * pa.GtK * pa.St * ( Vt - (fx.VtK(1,i)/Sa_p ) ) / pa.F;
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % -----> Intracellular Equations -----> 
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            
            fx.Jw = fx.Qb - sum(fx.Qa);
            fx.Nai = (1/var.w)*(fx.JNkcc1 - 3* fx.JNaK + fx.JNhe1 - fx.JAe4 - fx.Jw * var.Na);
            fx.Ki = (1/var.w)*(fx.JNkcc1 +2 * fx.JNaK - fx.JK - fx.Jw * var.K);
            fx.Cli = (1/var.w)*(2*fx.JNkcc1 + sum(fx.JCl) + fx.JAe4 - fx.Jw * var.Cl);
            fx.HCO3i = (1/var.w)*(fx.JBB - 2*fx.JAe4 - fx.Jw * var.HCO3);
            fx.Hyi = (1/var.w)*(fx.JBB - fx.JNhe1 - fx.Jw * var.H);
            
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
            % --------------- END -----------------
            %%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%
        end 
    end
end

