
%%%%%%%Setup value for simulink model%%%%%%

%%%%%%%%%%%%Generator Values%%%%%%%%%%%%%%%
%Calc base reactance
Pgen = 2500e3; %gen rating, don't confuse w. alternator rating
Vnom = 416; 
Zgen = Vnom.^2/Pgen;

%Reactances (pu)
Xd = 2.655; X1d = 0.186; X2d = 0.137; 
Xq = 1.766; X2q = 0.256; XL = 0.080;
Xnegseq = 0.197; Xzeroseq = 0.027;
%copypaste: [Xd X1d X2d Xq X2q XL]

%Resistances (converted from ohm to pu)
Rstatorwinding = 0.000256/Zgen;
Rrotorwinding = 1.630/Zgen;
Rexciterstatorfield = 17/Zgen;
Rexciterrotor = 0.092/Zgen;
Rpmgstator = 3.8/Zgen;

%Characteristic voltages
Vexcitationnoload = 15.0;
Vexcitationfullload = 67.0;

%Time constants
T1d = 0.213;
T2d = 0.016;
T1dOC = 5.1;
T2qo = T2d; %approximation used by Cummins
Tarmature = 0.081;
%copypaste: [T1d T2d T2qo]

%Physical parameters
Jgen = 125; %inertia, kg/m2
Hgen = 0.5*125*188.5^2/Pgen;  %inertia coefficient

%%%%%%%%%%%AVR Model parameters%%%%%%%%%%%%%
%Regulator constants
Tr = 0.018; Kpr = 0.998; Kir = 1; Kdr = 0.1; Tdr = 0.018;
Ka = 27.417; Ta = 0.019;
Kc = 0.1; Kd = 1; Ke = 1; Te = 0.3;

%PID/Regulator limits
Vpidmax = 10; Vpidmin = -10;
Vrmax = 10; Vrmin = 0;
Vfemax = 8; Vemin = 0;

%Saturation
E1 = 6; 
Se1 = 0.13;
E2 = 8;
Se2 = 0.82;

%AVR input filter zero/pole
Tfz = -10.612;
Tfp = -9.990;

%set up Gaussian filter for output filtering
%sigma = 5;
%size = 30;
%x = linspace(-size / 2, size / 2, size);
%gaussFilter = exp(-x .^ 2 / (3 * sigma ^ 2));
%gaussFilter = gaussFilter / sum (gaussFilter); % normalize

%%%%%%%%%%%%%%%Transformer Values%%%%%%%%%%%%%
%set up saturation values

%as given by manufacturer:
saturation_data = [0	0;
            0.0004	0.4008695332;
            0.0009	0.6013042998;
            0.0013	0.7015216831;
            0.0017	0.8017390664;
            0.0026	0.9019564497;
            0.0055	1.0021738329;
            0.0332	1.1023912162;
            0.3291	1.2026085995;];
        
%truncated around zero current for extrapolation function in Tx model
saturation_data_approx = [0 0; 
            0       0.7015;
            .0017	0.8017390664;
            0.0026	0.9019564497;
            0.0055	1.0021738329;
            0.0332	1.1023912162;
            0.3291	1.2026085995;];

        
%saturation_list is parameter in Tx blockset
saturation_list = saturation_data_approx;       
        
%initialise random initial fluxes:
%flux = 2.*rand(3,1)-1
%however, ABB says max. remanent flux in region 0.6-0.7 of sat.flux
%can't have 1p.u. remanent flux
%hence 0.7*(2.*rand(3,1)-1)
 
TxPwrRating = 2500e3; %VA
VTxHi = 24900; %V
VTxLo = 416; %V

%rest of parameters are pu
%RTx1 = 0.0048386; %original
RTx1 = 10.8/3/(VTxHi^2/TxPwrRating); %converted from datasheet
RTx2 = 0.0006117/3/(VTxLo^2/TxPwrRating);
%LTx2 = 0.0066061; %original
%RTx2 = 0.0006117;
LTx1 = 10.8*10.2/3/(VTxHi^2/TxPwrRating); %conv from datasheet
LTx2 = 10.2*RTx2;
MagRes = 500;
MagInd = 500;
%define remanent flux individually (if need be)
fluxlist1 = [-1 -1 1]; fluxlist2 = [-1 -1 1];
fluxlist3 = [-1 -1 1]; fluxlist4 = [-1 -1 1];
fluxlist5 = [-1 -1 1]; fluxlist6 = [-1 -1 1];
fluxlist7 = [-1 -1 1]; fluxlist8 = [-1 -1 1];
fluxlist9 = [-1 -1 1]; fluxlist10 = [-1 -1 1];
fluxlist11 = [-1 -1 1]; fluxlist12 = [-1 -1 1];
fluxlist13 = [-1 -1 1]; fluxlist14 = [-1 -1 1];
fluxlist15 = [-1 -1 1]; fluxlist16 = [-1 -1 1];
fluxlist17 = [-1 -1 1];

%%%%%%%%%%%%Cable parameters%%%%%%%%%%%
lcable = 250/3.28/1000; %cable length in km
rcable = (1.039*0.1438)/3.28; %resistance per km
hcable = 0.01880*3.28/(2*pi*60); %inductance per km
%fcable = %capacitance per km

%cposcable = 
%c0cable =

%%%%%Run model
%tic
%SimOut = sim('AVR_PLayground');
%toc

