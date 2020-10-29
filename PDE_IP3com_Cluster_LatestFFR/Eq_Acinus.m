function F = Eq_Acinus(c)
Nal_1  = c(1).var.Nal(1);
%%%%%
Nal_21 = c(2).var.Nal(1);
Nal_2  = c(2).var.Nal(2);
%%%%%
Nal_31 = c(3).var.Nal(1);
Nal_32 = c(3).var.Nal(2);
Nal_3  = c(3).var.Nal(3);
%%%%%
Nal_41 = c(4).var.Nal(1);
Nal_42 = c(4).var.Nal(2);
Nal_43 = c(4).var.Nal(3);
Nal_4  = c(4).var.Nal(4);
%%%%%
Nal_52 = c(5).var.Nal(1);
Nal_53 = c(5).var.Nal(2);
Nal_54 = c(5).var.Nal(3);
%%%%%
Nal_62 = c(6).var.Nal(1);
Nal_64 = c(6).var.Nal(2);
Nal_65 = c(6).var.Nal(3);
Nal_6  = c(6).var.Nal(4);
%%%%%%
Nal_75 = c(7).var.Nal(1);
Nal_76 = c(7).var.Nal(2);
%%%%%

Kl_1  = c(1).var.Kl(1);
%%%%%
Kl_21 = c(2).var.Kl(1);
Kl_2  = c(2).var.Kl(2);
%%%%%
Kl_31 = c(3).var.Kl(1);
Kl_32 = c(3).var.Kl(2);
Kl_3  = c(3).var.Kl(3);
%%%%%
Kl_41 = c(4).var.Kl(1);
Kl_42 = c(4).var.Kl(2);
Kl_43 = c(4).var.Kl(3);
Kl_4  = c(4).var.Kl(4);
%%%%%
Kl_52 = c(5).var.Kl(1);
Kl_53 = c(5).var.Kl(2);
Kl_54 = c(5).var.Kl(3);
%%%%%
Kl_62 = c(6).var.Kl(1);
Kl_64 = c(6).var.Kl(2);
Kl_65 = c(6).var.Kl(3);
Kl_6  = c(6).var.Kl(4);
%%%%%%
Kl_75 = c(7).var.Kl(1);
Kl_76 = c(7).var.Kl(2);
%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% First Calculate all the Cellular Concentrations.
%%%%% Cell 7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F(1)  = c(7).fx.Jw;
F(2)  = c(7).fx.Nai;
F(3)  = c(7).fx.Ki;
F(4)  = c(7).fx.Cli;
F(5)  = c(7).fx.HCO3i;
F(6)  = c(7).fx.Hyi;
%%%%% Cell 6
F(7)  = c(6).fx.Jw;
F(8)  = c(6).fx.Nai;
F(9)  = c(6).fx.Ki;
F(10) = c(6).fx.Cli;
F(11) = c(6).fx.HCO3i;
F(12) = c(6).fx.Hyi;
%%%%% Cell 5
F(13) = c(5).fx.Jw;
F(14) = c(5).fx.Nai;
F(15) = c(5).fx.Ki;
F(16) = c(5).fx.Cli;
F(17) = c(5).fx.HCO3i;
F(18) = c(5).fx.Hyi;
%%%%% Cell 4
F(19) = c(4).fx.Jw;
F(20) = c(4).fx.Nai;
F(21) = c(4).fx.Ki;
F(22) = c(4).fx.Cli;
F(23) = c(4).fx.HCO3i;
F(24) = c(4).fx.Hyi;
%%%%% Cell 3
F(25) = c(3).fx.Jw;
F(26) = c(3).fx.Nai;
F(27) = c(3).fx.Ki;
F(28) = c(3).fx.Cli;
F(29) = c(3).fx.HCO3i;
F(30) = c(3).fx.Hyi;
%%%%% Cell 2
F(31) = c(2).fx.Jw;
F(32) = c(2).fx.Nai;
F(33) = c(2).fx.Ki;
F(34) = c(2).fx.Cli;
F(35) = c(2).fx.HCO3i;
F(36) = c(2).fx.Hyi;
%%%%% Cell 1
F(37) = c(1).fx.Jw;
F(38) = c(1).fx.Nai;
F(39) = c(1).fx.Ki;
F(40) = c(1).fx.Cli;
F(41) = c(1).fx.HCO3i;
F(42) = c(1).fx.Hyi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%
% Lumen shared between cell 7 and cell 5.
%%%%%
JtNa_75 = c(7).fx.JtNa(1) + c(5).fx.JtNa(5);
JtK_75  = c(7).fx.JtK(1)  + c(5).fx.JtK(5);
QT_75   = c(7).fx.Qtot(1) + c(5).fx.Qtot(5);

F(43) = JtNa_75 - QT_75 * Nal_75;
F(44) = JtK_75  - QT_75 * Kl_75;

%%%%%
% Lumen shared between cell 7 and cell 6.
%%%%%

JtNa_76 = c(7).fx.JtNa(2) + c(6).fx.JtNa(5);
JtK_76  = c(7).fx.JtK(2)  + c(6).fx.JtK(5);
QT_76   = c(7).fx.Qtot(2) + c(6).fx.Qtot(5);

F(45) = JtNa_76 - QT_76 * Nal_76;
F(46) = JtK_76  - QT_76 * Kl_76;


%%%%%
% There is a bit of lumen which isn't shared and it comes from
% cell 6. This trickles down onto the first junction.
%%%%%
QT_6 = c(6).fx.Qtot(4);

F(47) =  c(6).fx.JtNa(4) - QT_6 * Nal_6;
F(48) =  c(6).fx.JtK(4)  - QT_6 * Kl_6;

%%%%%
% Lumen shared between cell 6 and 2.
%%%%%
JtNa_62 = c(6).fx.JtNa(1) + c(2).fx.JtNa(6);
JtK_62  = c(6).fx.JtK(1)  + c(2).fx.JtK(6);
QT_62   = c(6).fx.Qtot(1) + c(2).fx.Qtot(6);

F(49) = JtNa_62 - QT_62 * Nal_62;
F(50) = JtK_62  - QT_62 * Kl_62;

%%%%%
% First Junction. Space shared between cell 6 and 5.
%%%%%

JtNa_65 = c(6).fx.JtNa(3) + c(5).fx.JtNa(4);
JtK_65  = c(6).fx.JtK(3)  + c(5).fx.JtK(4);
QT_65 = QT_75 + QT_76 + QT_6 + QT_62 + c(6).fx.Qtot(3) + c(5).fx.Qtot(4);


F(51) = JtNa_65 ...
        + QT_75 * Nal_75 ...
        + QT_76 * Nal_76 ...
        + QT_62 * Nal_62 ...
        + QT_6  * Nal_6 ...
        - QT_65 * Nal_65;    

F(52) = JtK_65 ...
        + QT_75 * Kl_75 ...
        + QT_76 * Kl_76 ...
        + QT_62 * Kl_62 ...
        + QT_6  * Kl_6 ...
        - QT_65 * Kl_65;

%%%%%%
% Luminal Space shared between cell 2 and 5.
%%%%%%

JtNa_52 = c(5).fx.JtNa(1) + c(2).fx.JtNa(5);
JtK_52  = c(5).fx.JtK(1)  + c(2).fx.JtK(5);
QT_52   = QT_65 + c(5).fx.Qtot(1) + c(2).fx.Qtot(5);

F(53) = JtNa_52 ...
    + QT_65 * Nal_65 ...
    - QT_52 * Nal_52;

F(54) = JtK_52 ...
    + QT_65 * Kl_65 ...
    - QT_52 * Kl_52;

%%%%%%
% Luminal Space shared between cell 5 and 4.
%%%%%%

JtNa_54 = c(5).fx.JtNa(3)+ c(4).fx.JtNa(5);
JtK_54  = c(5).fx.JtK(3) + c(4).fx.JtK(5);
QT_54   = c(5).fx.Qtot(3)+ c(4).fx.Qtot(5);

F(55) = JtNa_54 ...
    - QT_54 * Nal_54;

F(56) = JtK_54 ...
    - QT_54 * Kl_54;

%%%%%%%%
% Luminal Space shared between cell 6 and 4.
%%%%%%%%

JtNa_64 = c(6).fx.JtNa(2)+ c(4).fx.JtNa(6);
JtK_64  = c(6).fx.JtK(2) + c(4).fx.JtK(6);
QT_64   = c(6).fx.Qtot(2)+ c(4).fx.Qtot(6);

F(57) = JtNa_64 - QT_64 * Nal_64;
F(58) = JtK_64  - QT_64 * Kl_64;

%%%%%%%%
% Luminal Space shared between cell 4 and 2.
%%%%%%%%

JtNa_42 = c(4).fx.JtNa(2) + c(2).fx.JtNa(4);
JtK_42 = c(4).fx.JtK(2) + c(2).fx.JtK(4);
QT_42 = QT_64 + c(4).fx.Qtot(2) + c(2).fx.Qtot(4);

F(59) = JtNa_42 ...
    + QT_64 * Nal_64 ...
    - QT_42 * Nal_42;

F(60) = JtK_42 ...
    + QT_64 * Kl_64 ...
    - QT_42 * Kl_42;


%%%%%%%%
% Luminal Space of cell 4 .
%%%%%%%%

JtNa_4 = c(4).fx.JtNa(4);
JtK_4 = c(4).fx.JtK(4);
QT_4 = c(4).fx.Qtot(4);

F(61) = JtNa_4 ...
    - QT_4 * Nal_4;

F(62) = JtK_4 ...
    - QT_4 * Kl_4;


%%%%%%%%
% Luminal Space shared between cell 4 and 3.
%%%%%%%%

JtNa_43 = c(4).fx.JtNa(3) + c(3).fx.JtNa(4);
JtK_43 = c(4).fx.JtK(3) + c(3).fx.JtK(4);
QT_43 = QT_4 + c(4).fx.Qtot(3) + c(3).fx.Qtot(4);

F(63) = JtNa_43 ...
    + QT_4 * Nal_4 ...
    - QT_43 * Nal_43;

F(64) = JtK_43 ...
    + QT_4 * Kl_4 ...
    - QT_43 * Kl_43;


%%%%%%%%
% Luminal Space of Cell 1 .
%%%%%%%%

JtNa_1 = c(1).fx.JtNa(1);
JtK_1  = c(1).fx.JtK(1);
QT_1   = c(1).fx.Qtot(1);

F(65) = JtNa_1 ...
    - QT_1 * Nal_1;

F(66) = JtK_1 ...
    - QT_1 * Kl_1;

%%%%%%%%
% Luminal Space shared by Cells 3 and 1 .
%%%%%%%%

JtNa_31 = c(3).fx.JtNa(1) + c(1).fx.JtNa(3);
JtK_31  = c(3).fx.JtK(1)  + c(1).fx.JtK(3);
QT_31   = QT_1 + c(3).fx.Qtot(1) + c(1).fx.Qtot(3);

F(67) = JtNa_31 ...
    + QT_1 * Nal_1 ...
    - QT_31 * Nal_31;

F(68) = JtK_31 ...
    + QT_1 * Kl_1 ...
    - QT_31 * Kl_31;

%%%%%%%%
% Luminal Space shared by Cells 4 and 1 .
%%%%%%%%

JtNa_41 = c(4).fx.JtNa(1) + c(1).fx.JtNa(4);
JtK_41  = c(4).fx.JtK(1)  + c(1).fx.JtK(4);
QT_41   = c(4).fx.Qtot(1) + c(1).fx.Qtot(4);

F(69) = JtNa_41 ...
    - QT_41 * Nal_41;

F(70) = JtK_41 ...
    - QT_41 * Kl_41;

%%%%%%%%
% Luminal Space shared by Cells 2 and 1 .
%%%%%%%%

JtNa_21 = c(2).fx.JtNa(1) + c(1).fx.JtNa(2);
JtK_21  = c(2).fx.JtK(1)  + c(1).fx.JtK(2);
QT_21   = QT_41 + QT_31 + c(2).fx.Qtot(1) + c(1).fx.Qtot(2);

F(71) = JtNa_21 ...
    + QT_41 * Nal_41 ...
    + QT_31 * Nal_31 ...
    - QT_21 * Nal_21;

F(72) = JtK_21 ...
    + QT_41 * Kl_41 ...
    + QT_31 * Kl_31 ...
    - QT_21 * Kl_21;  

%%%%%%%%
% Luminal Space of Cell 2 .
%%%%%%%%

JtNa_2 = c(2).fx.JtNa(2);
JtK_2  = c(2).fx.JtK(2);
QT_2   = c(2).fx.Qtot(2);

F(73) = JtNa_2 ...
    - QT_2 * Nal_2;

F(74) = JtK_2 ...
    - QT_2 * Kl_2;

%%%%%%%%
% Luminal Space shared by Cells 5 and 3 .
%%%%%%%%

JtNa_53 = c(5).fx.JtNa(2) + c(3).fx.JtNa(5);
JtK_53  = c(5).fx.JtK(2)  + c(3).fx.JtK(5);
QT_53   = c(5).fx.Qtot(2) + c(3).fx.Qtot(5);

F(75) = JtNa_53 ...
    - QT_53 * Nal_53;

F(76) = JtK_53 ...
    - QT_53 * Kl_53;

%%%%%%%%
% Luminal Space shared by Cells 5 and 3 .
%%%%%%%%

JtNa_32 = c(3).fx.JtNa(2) + c(2).fx.JtNa(3);
JtK_32  = c(3).fx.JtK(2)  + c(2).fx.JtK(3);
QT_32   = QT_54 + QT_42 + QT_43 + QT_21 ...
        + QT_2 + QT_53 ...
        + c(3).fx.Qtot(2) + c(2).fx.Qtot(3);

F(77) = JtNa_32 ...
    + QT_54 * Nal_54 ...
    + QT_42 * Nal_42 ...
    + QT_43 * Nal_43 ...
    + QT_21 * Nal_21 ...
    + QT_2  * Nal_2  ...
    + QT_53 * Nal_53 ...
    - QT_32 * Nal_32;

F(78) = JtK_32 ...
    + QT_54 * Kl_54 ...
    + QT_42 * Kl_42 ...
    + QT_43 * Kl_43 ...
    + QT_21 * Kl_21 ...
    + QT_2  * Kl_2  ...
    + QT_53 * Kl_53 ...
    - QT_32 * Kl_32;

%%%%%%%%
% Luminal Space of Cell 3. EXIT
%%%%%%%%

JtNa_3 = c(3).fx.JtNa(3);
JtK_3  = c(3).fx.JtK(3);
QT_3   = QT_32 + c(3).fx.Qtot(3);

F(79) = JtNa_3 ...
    + QT_32 * Nal_32 ...
    - QT_3 * Nal_3;

F(80) = JtK_3 ...
    + QT_32 * Kl_32 ...
    - QT_3 * Kl_3;

F = F';
end