w1    = U(:,1);
Na1   = U(:,2);
K1    = U(:,3);
Cl1   = U(:,4);
HCO31 = U(:,5);
H1    = U(:,6);
Va1   = U(:,7);
Vb1   = U(:,8);
w2    = U(:,9);
Na2   = U(:,10);
K2    = U(:,11);
Cl2   = U(:,12);
HCO32 = U(:,13);
H2    = U(:,14);
Va2   = U(:,15);
Vb2   = U(:,16);
w3    = U(:,17);
Na3   = U(:,18);
K3    = U(:,19);
Cl3   = U(:,20);
HCO33 = U(:,21);
H3    = U(:,22);
Va3   = U(:,23);
Vb3   = U(:,24);
w4    = U(:,25);
Na4   = U(:,26);
K4    = U(:,27);
Cl4   = U(:,28);
HCO34 = U(:,29);
H4    = U(:,30);
Va4   = U(:,31);
Vb4   = U(:,32);
w5    = U(:,33);
Na5   = U(:,34);
K5    = U(:,35);
Cl5   = U(:,36);
HCO35 = U(:,37);
H5    = U(:,38);
Va5   = U(:,39);
Vb5   = U(:,40);
w6    = U(:,41);
Na6   = U(:,42);
K6    = U(:,43);
Cl6   = U(:,44);
HCO36 = U(:,45);
H6    = U(:,46);
Va6   = U(:,47);
Vb6   = U(:,48);
w7    = U(:,49);
Na7   = U(:,50);
K7    = U(:,51);
Cl7   = U(:,52);
HCO37 = U(:,53);
H7    = U(:,54);
Va7   = U(:,55);
Vb7   = U(:,56);



Nal11 = U(:,57);
Nal21 = U(:,58);
Nal22 = U(:,59);
Nal31 = U(:,60);
Nal32 = U(:,61);
Nal33 = U(:,62);
Nal41 = U(:,63);
Nal42 = U(:,64);
Nal43 = U(:,65);
Nal44 = U(:,66);
Nal52 = U(:,67);
Nal53 = U(:,68);
Nal54 = U(:,69);
Nal62 = U(:,70);
Nal64 = U(:,71);
Nal65 = U(:,72);
Nal66 = U(:,73);
Nal75 = U(:,74);
Nal76 = U(:,75);
Kl75  = U(:,76);
Kl76  = U(:,77);

Kl11  = U(:,78);
Kl21  = U(:,79);
Kl22  = U(:,80);
Kl31  = U(:,81);
Kl32  = U(:,82);
Kl33  = U(:,83);
Kl41  = U(:,84);
Kl42  = U(:,85);
Kl43  = U(:,86);
Kl44  = U(:,87);
Kl52  = U(:,88);
Kl53  = U(:,89);
Kl54  = U(:,90);
Kl62  = U(:,91);
Kl64  = U(:,92);
Kl65  = U(:,93);
Kl66  = U(:,94);

Cll11 = U(:,95);
Cll21 = U(:,96);
Cll22 = U(:,97);
Cll31 = U(:,98);
Cll32 = U(:,99);
Cll33 = U(:,100);
Cll41 = U(:,101);
Cll42 = U(:,102);
Cll43 = U(:,103);
Cll44 = U(:,104);
Cll52 = U(:,105);
Cll53 = U(:,106);
Cll54 = U(:,107);
Cll62 = U(:,108);
Cll64 = U(:,109);
Cll65 = U(:,110);
Cll66 = U(:,111);
Cll75 = U(:,112);
Cll76 = U(:,113);




Ii = 2*(Na1 + K1 + H1) + par.CO20;

Il11 = 2*Cll11 + par.Ul;
Il12 = 2*Cll21 + par.Ul;
Il13 = 2*Cll31 + par.Ul;
Il14 = 2*Cll41 + par.Ul;

Qa11 = par.Sa_p{1,1} * par.La * (Il11 - Ii);
Qa12 = par.Sa_p{1,2} * par.La * (Il12 - Ii);
Qa13 = par.Sa_p{1,3} * par.La * (Il13 - Ii);
Qa14 = par.Sa_p{1,4} * par.La * (Il14 - Ii);

Qt11 = par.Sa_p{1,1} * par.Lt * (Il11 - par.Ie);
Qt12 = par.Sa_p{1,2} * par.Lt * (Il12 - par.Ie);
Qt13 = par.Sa_p{1,3} * par.Lt * (Il13 - par.Ie);
Qt14 = par.Sa_p{1,4} * par.Lt * (Il14 - par.Ie);

Qtot11 =(Qa11 + Qt11);
Qtot12 =(Qa12 + Qt12);
Qtot13 =(Qa13 + Qt13);
Qtot14 =(Qa14 + Qt14);


Ii2 = 2*(Na2 + K2 + H2) + par.CO20;
    
Il21 = 2*Cll21 + par.Ul;
Il22 = 2*Cll22 + par.Ul;
Il23 = 2*Cll32 + par.Ul;
Il24 = 2*Cll42 + par.Ul;
Il25 = 2*Cll52 + par.Ul;
Il26 = 2*Cll62 + par.Ul;

Qa21 = par.Sa_p{2,1} * par.La * (Il21 - Ii2);
Qa22 = par.Sa_p{2,2} * par.La * (Il22 - Ii2);
Qa23 = par.Sa_p{2,3} * par.La * (Il23 - Ii2);
Qa24 = par.Sa_p{2,4} * par.La * (Il24 - Ii2);
Qa25 = par.Sa_p{2,5} * par.La * (Il25 - Ii2);
Qa26 = par.Sa_p{2,6} * par.La * (Il26 - Ii2);

Qt21 = par.Sa_p{2,1} * par.Lt * ( Il21 - par.Ie );
Qt22 = par.Sa_p{2,2} * par.Lt * ( Il22 - par.Ie );
Qt23 = par.Sa_p{2,3} * par.Lt * ( Il23 - par.Ie );
Qt24 = par.Sa_p{2,4} * par.Lt * ( Il24 - par.Ie );
Qt25 = par.Sa_p{2,5} * par.Lt * ( Il25 - par.Ie );
Qt26 = par.Sa_p{2,6} * par.Lt * ( Il26 - par.Ie );

Qtot21 =(Qa21 + Qt21);
Qtot22 =(Qa22 + Qt22);
Qtot23 =(Qa23 + Qt23);
Qtot24 =(Qa24 + Qt24);
Qtot25 =(Qa25 + Qt25);
Qtot26 =(Qa26 + Qt26);

Ii3 = 2*(Na3 + K3 + H3) + par.CO20;
    
Il31 = 2*Cll31 + par.Ul;
Il32 = 2*Cll32 + par.Ul;
Il33 = 2*Cll33 + par.Ul;
Il34 = 2*Cll43 + par.Ul;
Il35 = 2*Cll53 + par.Ul;

Qa31 = par.Sa_p{3,1} * par.La * (Il31 - Ii3);
Qa32 = par.Sa_p{3,2} * par.La * (Il32 - Ii3);
Qa33 = par.Sa_p{3,3} * par.La * (Il33 - Ii3);
Qa34 = par.Sa_p{3,4} * par.La * (Il34 - Ii3);
Qa35 = par.Sa_p{3,5} * par.La * (Il35 - Ii3);

Qt31 = par.Sa_p{3,1} * par.Lt * (Il31 - par.Ie);
Qt32 = par.Sa_p{3,2} * par.Lt * (Il32 - par.Ie);
Qt33 = par.Sa_p{3,3} * par.Lt * (Il33 - par.Ie);
Qt34 = par.Sa_p{3,4} * par.Lt * (Il34 - par.Ie);
Qt35 = par.Sa_p{3,5} * par.Lt * (Il35 - par.Ie);

Qtot31 =(Qa31 + Qt31);
Qtot32 =(Qa32 + Qt32);
Qtot33 =(Qa33 + Qt33);
Qtot34 =(Qa34 + Qt34);
Qtot35 =(Qa35 + Qt35);
   

Ii4 = 2*(Na4 + K4 + H4) + par.CO20;

Il41 = 2*Cll41 + par.Ul;
Il42 = 2*Cll42 + par.Ul;
Il43 = 2*Cll43 + par.Ul;
Il44 = 2*Cll44 + par.Ul;
Il45 = 2*Cll54 + par.Ul;
Il46 = 2*Cll64 + par.Ul;

Qa41 = par.Sa_p{4,1} * par.La * (Il41 - Ii4);
Qa42 = par.Sa_p{4,2} * par.La * (Il42 - Ii4);
Qa43 = par.Sa_p{4,3} * par.La * (Il43 - Ii4);
Qa44 = par.Sa_p{4,4} * par.La * (Il44 - Ii4);
Qa45 = par.Sa_p{4,5} * par.La * (Il45 - Ii4);
Qa46 = par.Sa_p{4,6} * par.La * (Il46 - Ii4);

Qt41 = par.Sa_p{4,1} * par.Lt * (Il41 - par.Ie);
Qt42 = par.Sa_p{4,2} * par.Lt * (Il42 - par.Ie);
Qt43 = par.Sa_p{4,3} * par.Lt * (Il43 - par.Ie);
Qt44 = par.Sa_p{4,4} * par.Lt * (Il44 - par.Ie);
Qt45 = par.Sa_p{4,5} * par.Lt * (Il45 - par.Ie);
Qt46 = par.Sa_p{4,6} * par.Lt * (Il46 - par.Ie);

Qtot41 =(Qa41 + Qt41);
Qtot42 =(Qa42 + Qt42);
Qtot43 =(Qa43 + Qt43);
Qtot44 =(Qa44 + Qt44);
Qtot45 =(Qa45 + Qt45);
Qtot46 =(Qa46 + Qt46);


Ii5 = 2*(Na5 + K5 + H5) + par.CO20;
    
Il52 = 2*Cll52 + par.Ul;
Il53 = 2*Cll53 + par.Ul;
Il54 = 2*Cll54 + par.Ul;
Il56 = 2*Cll65 + par.Ul;
Il57 = 2*Cll75 + par.Ul;

Qa52 = par.Sa_p{5,2} * par.La * (Il52 - Ii5);
Qa53 = par.Sa_p{5,3} * par.La * (Il53 - Ii5);
Qa54 = par.Sa_p{5,4} * par.La * (Il54 - Ii5);
Qa56 = par.Sa_p{5,6} * par.La * (Il56 - Ii5);
Qa57 = par.Sa_p{5,7} * par.La * (Il57 - Ii5);

Qt52 = par.Sa_p{5,2} * par.Lt * ( Il52 - par.Ie );
Qt53 = par.Sa_p{5,3} * par.Lt * ( Il53 - par.Ie );
Qt54 = par.Sa_p{5,4} * par.Lt * ( Il54 - par.Ie );
Qt56 = par.Sa_p{5,6} * par.Lt * ( Il56 - par.Ie );
Qt57 = par.Sa_p{5,7} * par.Lt * ( Il57 - par.Ie );

Qtot52 =(Qa52 + Qt52);
Qtot53 =(Qa53 + Qt53);
Qtot54 =(Qa54 + Qt54);
Qtot56 =(Qa56 + Qt56);
Qtot57 =(Qa57 + Qt57);


Ii6 = 2*(Na6 + K6 + H6) + par.CO20;
Il62 = 2*Cll62 + par.Ul;
Il64 = 2*Cll64 + par.Ul;
Il65 = 2*Cll65 + par.Ul;
Il66 = 2*Cll66 + par.Ul;
Il67 = 2*Cll76 + par.Ul;

Qa62 = par.Sa_p{6,2} * par.La * (Il62 - Ii6);
Qa64 = par.Sa_p{6,4} * par.La * (Il64 - Ii6);
Qa65 = par.Sa_p{6,5} * par.La * (Il65 - Ii6);
Qa66 = par.Sa_p{6,6} * par.La * (Il66 - Ii6);
Qa67 = par.Sa_p{6,7} * par.La * (Il67 - Ii6);

Qt62 = par.Sa_p{6,2} * par.Lt * ( Il62 - par.Ie );
Qt64 = par.Sa_p{6,4} * par.Lt * ( Il64 - par.Ie );
Qt65 = par.Sa_p{6,5} * par.Lt * ( Il65 - par.Ie );
Qt66 = par.Sa_p{6,6} * par.Lt * ( Il66 - par.Ie );
Qt67 = par.Sa_p{6,7} * par.Lt * ( Il67 - par.Ie );

Qtot62 =(Qa62 + Qt62);
Qtot64 =(Qa64 + Qt64);
Qtot65 =(Qa65 + Qt65);
Qtot66 =(Qa66 + Qt66);
Qtot67 =(Qa67 + Qt67);

Ii7 = 2*(Na7 + K7 + H7) + par.CO20;
Il75 = 2*Cll75 + par.Ul;
Il76 = 2*Cll76 + par.Ul;

Qa75 = par.Sa_p{7,5} * par.La * (Il75 - Ii7);
Qa76 = par.Sa_p{7,6} * par.La * (Il76 - Ii7);
Qt75 = par.Sa_p{7,5} * par.Lt * ( Il75 - par.Ie );
Qt76 = par.Sa_p{7,6} * par.Lt * ( Il76 - par.Ie ); 

Qtot75 =(Qa75 + Qt75);
Qtot76 =(Qa76 + Qt76);

Q75 = (Qtot75 + Qtot57);
Q76 = (Qtot76 + Qtot67);
Q65 = (Qtot65 + Qtot56) + Q75 + Q76;
Q53 = (Qtot53 + Qtot35);
Q52 = (Qtot52 + Qtot25) + Q65 + Q53;
Q22 =  Qtot22;
Q32 = (Qtot32 + Qtot23) + Q22 + Q52;
Q62 = (Qtot62 + Qtot26);
Q64 = (Qtot64 + Qtot46) + Qtot66 + Qtot44;
Q42 = (Qtot42 + Qtot24) + Q62 + Q64;
Q54 = (Qtot54 + Qtot45);
Q43 = (Qtot43 + Qtot34) + Q54;
Q11 = Qtot11;
Q41 = (Qtot41 + Qtot14);
Q31 = (Qtot31 + Qtot13);
Q21 = (Qtot21 + Qtot12) + Q41 + Q31 + Q11;
Q33 =  Qtot33 + Q32 + Q43 + Q42 + Q21;
Q33(end)