function dx=Secretion(~,x,cav_tri,param, surf_tri,i,w_IPR,w_basal)

aNaK = param.aNaK(i);
aNkcc1 = param.aNkcc1(i);
GtNa = param.GtNa(i);
GtK  = param.GtK(i);
GCl = param.GCl(i);
GK = param.GK(i);
G1 = param.G1(i);
G4 = param.G4(i);
GB = param.GB(i);
St = param.St(i);
Sb = param.Sb(i);
Sa = param.Sa(i);

Nal = x(1);
Kl = x(2);
Cll = x(3);
w = x(4);
Na = x(5);
K = x(6);
Cl = x(7);
HCO3 = x(8);
H = x(9);
Va = x(10);
Vb = x(11);

%% Currents and fluxes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (3Na+)/(2K+) ATP-ase pump (NaK)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NaKbasalfactor = 0.7;                                                       % Fraction of NaK ATPase in the basal membrane

JNaKb = NaKbasalfactor*Sb*aNaK*( param.r*param.Ke^2 * Na^3 ... 
  / ( param.Ke^2 + param.alpha1 * Na^3 ) );   

JNaKa = (1-NaKbasalfactor)*Sa*aNaK*(param.r * Kl^2 * Na^3 ...
  / ( Kl^2 + param.alpha1 * Na^3 ) ); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nernst potentials (J/C)1e3 = mV
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

VCl = param.RTF * log( Cll / Cl );        
VKb = param.RTF * log( param.Ke / K );     
VKa = param.RTF * log( Kl / K );          
VtNa = param.RTF * log( Nal / param.Nae );
VtK = param.RTF * log( Kl / param.Ke );   
Vt = Va - Vb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ca2+ activated channels.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cl channels
PCl = 1 ./ ( 1 + ( param.KCaCC ./ (cav_tri.*(w_IPR>0)) ).^param.eta1 ); 
PrCl = sum( PCl .* surf_tri )/Sa;
JCl = GCl * PrCl * ( Va + VCl ) / param.F;          % fS.micro-metres^2.mV.mol.C^-1

% K Channels
PKb = 1 ./ ( 1 + ( param.KCaKC ./ (cav_tri.*(w_basal>0)) ).^param.eta2 ); 
PKa = 1 ./ ( 1 + ( param.KCaKC ./ (cav_tri.*(w_IPR>0)) ).^param.eta2 ); 
PrKb = sum( PKb .* surf_tri )/Sb;
PrKa = sum( PKa .* surf_tri )/Sa;


JKa = GK * PrKa * ( Va - VKa ) / param.F;            % fS.micro-metres^2.mV.mol.C^-1
JKb = GK * PrKb * ( Vb - VKb ) / param.F;            % fS.micro-metres^2.mV.mol.C^-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tight Junction Na+ and K+ currents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JtNa = GtNa * St * ( Vt - VtNa ) / param.F;                 % fS.micro-metres^2.mV.mol.C^-1
JtK = GtK * St * ( Vt - VtK ) / param.F;                      % fS.micro-metres^2.mV.mol.C^-1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Osmolarities 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qa = param.B1 * ( 2 * ( Nal + Kl - Na - K - H ) - param.CO20 + param.Ul );     % micro-metres^3.s^-1
Qb = param.B2 * ( 2 * ( Na + K + H ) + param.CO20 - ...
                      ( param.Nae + param.Ke + param.Cle + param.HCO3e ) );
Qt = param.B3 * ( 2 * ( Nal + Kl ) + param.Ul - ....
                      ( param.Nae + param.Ke + param.Cle + param.HCO3e ) ); % micro-metres^3.s^-1
Qtot=(Qa+Qt);                                     % micro-metres^3.s^-1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Na+ K+ 2Cl- co-transporter (Nkcc1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JNkcc1 = aNkcc1 * Sb * ( param.a1 - param.a2 * Na * K * Cl^2 ) ...
                                             / ( param.a3 + param.a4 * Na * K * Cl^2 );     
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (Na+)2 HCO3-/Cl- Anion exchanger (Ae4)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JAe4 = Sb * G4 * ( ( param.Cle / ( param.Cle + param.KCl ) ) * ( Na / ( Na + param.KNa ) ) ...
             * ( HCO3 / ( HCO3 + param.KB ) )^2 );       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Na+ / H+ Anion exchanger (Nhe1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JNhe1 = Sb * G1 * ( ( param.Nae / ( param.Nae + param.KNa ) ) * ( H / ( param.KH + H ) )...
                          - ( Na / ( Na + param.KNa ) ) * ( param.He / ( param.KH + param.He ) ) ); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bicarbonate Buffer (Reaction Term)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This equation is a reaction inside the cell, note how it depends on the
% cellular volume

JBB = w * GB * ( param.kp * param.CO20 - param.kn * HCO3 * H );                   
                  
%% Equations

dx(1) = ( JtNa - Qtot*Nal + 3*JNaKa )/param.wl(i);
dx(2) = ( JtK  - Qtot*Kl + JKa - 2*JNaKa)/param.wl(i);
dx(3) = ( -JCl - Qtot*Cll     )/param.wl(i);
dx(4) = Qb - Qa;
dx(5) = ( JNkcc1 - 3*(JNaKb+JNaKa) + JNhe1 - JAe4 - dx(4) * Na ) / w;
dx(6) = ( JNkcc1 + 2*(JNaKb+JNaKa) - JKb - JKa - dx(4) * K ) / w;
dx(7) = ( 2 * JNkcc1 + JAe4 + JCl - dx(4) * Cl ) / w;
dx(8) = ( JBB - 2 * JAe4 - dx(4) * HCO3 ) / w;
dx(9) = ( JBB - JNhe1 - dx(4) * H ) / w;
dx(10) = 100*(-JCl - JNaKa - JKa - JtK - JtNa);      % Note the arbitrary factor of 100, just to make sure Va is fast.
dx(11) = 100*(     - JNaKb - JKb + JtK + JtNa);

dx = dx';
end