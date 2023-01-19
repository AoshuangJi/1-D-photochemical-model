%% Section 1: FIGURE for OZONE COMPARISON AMONG THREE MODELS
%% 1.1 Ozone columd depth in the 1-D, WACCM-3D, AND ROCKE-3D MODELS
% ********** read data: Ozone column depth (DU) at different pO2 *********
pO2 = [1 0.1 0.01 0.001 0.0001]; % O2 levels

% ozone column depth in the 1d
% SZA = 48.2, AGL = 0.375, with old H2O cross section 
old_O3 = [268.15 210 107.35 18.098 1.9058]; 
% SZA = 48.2, AGL = 0.375, with new H2O cross section
new_O3 = [268.15 224.06 127.09 24.507 2.787]; 
% 1D gaussian calculation with high H2O at the cold trap
Gau_O3 = [287 236 128 24 3]; 
% Current standard 1-D model (8-pt, new H2O cross section, low H2O at the cold trap)
Gau_LOWH2O_O3 = [298 248 142 33 5];
% without lightning
O3_nolit=[250 201 126 24 3];
% 1D model used by Segura et al. in 2003
segura_O3 = [311 266 124 25.3 1.21];  
% evenly 8-pt
Even_O3 = [278 223 120 26];

%*** 1D model with chlorine***
% Standard run with CH3CL, fixed flux (2.92e8 cm-3 s-1) for low O2
O3_CL = [296 236 120 24]; 
% Low CH3CL (30 times lower CH3CL flux than that in the standard run)
O3_LOWCL = [297 249 141 33];
% Fixed CH3CL (0.5ppt) for all the pO2
O3_CL05ppt = [296 239 113 5];
O3_CL05pptNOSNOA = [306 263 13];

phi_CL = [2.92 2.12 4.40 28.8]*1E8;
phi_CLNOSNOA = [2.80 2.09 162 445]*1E8;
% Mixing ratios of CH3CL for standard 1d
mrCH3CL = [5.00E-10 6.81E-10 3.62E-10 1.26E-10 6.74E-11];

% ozone column depth in the WACCM 3D
Cooke_O3 = [280 165 65 15 ]; % mean
Cooke_O3_min = [235 130 40 10]; % minimum
Cooke_O3_max = [375 235 85 25]; % maxmum

% ozone column depth in the ROCKE 3D
Way_O3 = [298 247 123 30 2.5]; 

% ozone column depth in the Kinetics 1D model
Kinetics_O3 = [330 163 51 9]; 
%% ****** plot: O3 column depth vs. pO2 with 8 GauPt,Segura,&Cooke ******
%
figure
subplot(211)
plot(pO2,Gau_LOWH2O_O3,'o-','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(pO2(1:4),O3_CL,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(pO2,segura_O3,'o-','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(pO2(1:4),Kinetics_O3,'Color',[0 0 0],'LineWidth',2)
plot(pO2(1:4),Cooke_O3,'o-','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(pO2,Way_O3,'o-','Color',[0.5,0.3,0.8],'LineWidth',2)
plot(pO2,ones(5)*292,':','LineWidth',2)
scatter(pO2(1),375,70,'o','filled')
scatter(pO2(1),235,70,'o','filled')
scatter(pO2(2),235,70,'o','filled')
scatter(pO2(2),130,70,'o','filled')
scatter(pO2(3),85,70,'o','filled')
scatter(pO2(3),40,70,'o','filled')
scatter(pO2(4),25,70,'o','filled')
scatter(pO2(4),10,70,'o','filled')
errorbar(pO2(1:4),(Cooke_O3_min+Cooke_O3_max)/2,(Cooke_O3_max-Cooke_O3_min)/2)
text(1e-3,340,'WOUDC Average = 292 DU','fontsize',14)
ylim([0 400]);xlim([1e-4 1.5]);
yticks([0 50 100 150 200 250 300 350 400]);
xlabel('pO_2 (PAL)','fontweight','bold');ylabel('O_3 column depth (DU)','fontweight','bold');
legend('Standard 1D','Standard 1D with chlorine','Segura et al., 2003','Standard Kinetics 1D model','Cooke et al., 2022','Way et al., 2017','fontsize',14,'Location','best');
set(gca,'XScale','log');
set(gca,'LineWidth',2,'FontSize',20);
box on

subplot(212)
plot(pO2,Gau_LOWH2O_O3,'o-','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(pO2(1:4),O3_CL,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(pO2,segura_O3,'o-','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(pO2(1:4),Kinetics_O3,'Color',[0 0 0],'LineWidth',2)
plot(pO2(1:4),Cooke_O3,'o-','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(pO2,Way_O3,'o-','Color',[0.5,0.3,0.8],'LineWidth',2)
plot(pO2,ones(5)*292,':','LineWidth',2)
xlim([1e-4 1.5]);
text(1e-3,340,'WOUDC Average = 292 DU','fontsize',14)
xlabel('pO_2 (PAL)','fontweight','bold');ylabel('O_3 column depth (DU)','fontweight','bold');
set(gca,'XScale','log');set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
box on

%% Read data from 1D outputs
% data includes
% Z:Altitude in cm. To get km:./1e5
% PRESS: Pressure. To get hPa or mbar: ./1e3
% DEN: Total number density in cm^-3.
% ProOH1: Production rate of OH from O1D + H2O in cm^-3 s^-1.
% ProOH25: Production rate of OH from H2O + HV in cm^-3 s^-1.
% ProO3D: Production rate of O1D from O3 + HV -> O2 + O(1D) in cm^-3 s^-1.
% PO2D: Photolysis rate of O2 + HV -> O + O(1D) in s^-1.
% PO2: Photolysis rate of O2 + HV -> O + O(3P) in s^-1.
% PH2O: Photolysis rate of H2O + HV -> H + OH in s^-1.
% PO3D: Photolysis rate of O3 + HV -> O2 + O(1D) in s^-1.
% PH2O2: Photolysis rate of H2O2 + HV -> H + HO2 in s^-1.
% FO3: Mixing ratio of O3
% FH2O: Mixing ratio of H2O
% NumO2: Number density of O2 in cm^-3.
% NumOH: Number density of OH in cm^-3.
% NumO1D: Number density of O1D in cm^-3.

% find the pathway
dir1 = '/Users/aoshuang/Documents/MATLAB/PHOTOMODEL/OUTPUT_PLOT0828/';
dir2 = 'Converged/plot/';
dir3 = 'ONESTEP/plot/';

%Converged Data for standard
file1 = fullfile(dir1,dir2,'OUTPUT_PLOT1PAL.dat');
file2 = fullfile(dir1,dir2,'OUTPUT_PLOT01PAL.dat');
file3 = fullfile(dir1,dir2,'OUTPUT_PLOT001PAL.dat');
file4 = fullfile(dir1,dir2,'OUTPUT_PLOT0001PAL.dat');

DATA_1PAL_1D = readmatrix(file1);
DATA_01PAL_1D = readmatrix(file2);
DATA_001PAL_1D = readmatrix(file3);
DATA_0001PAL_1D = readmatrix(file4);
% Converged Data without scattering 
file1_NOS = fullfile(dir1,dir2,'OUTPUT_PLOT1PAL_NOS.dat');
file2_NOS = fullfile(dir1,dir2,'OUTPUT_PLOT01PAL_NOS.dat');
file3_NOS = fullfile(dir1,dir2,'OUTPUT_PLOT001PAL_NOS.dat');
file4_NOS = fullfile(dir1,dir2,'OUTPUT_PLOT0001PAL_NOS.dat');
file0920_1 = fullfile(dir1,'OUTPUT_PLOT0001PAL_NOS0920.dat');
file0920_2 = fullfile(dir1,'OUTPUT_PLOT0001PAL_NOS_SR.dat');

Data_1PAL_1D_NOS = readmatrix(file1_NOS);
Data_01PAL_1D_NOS = readmatrix(file2_NOS);
Data_001PAL_1D_NOS = readmatrix(file3_NOS);
Data_0001PAL_1D_NOS = readmatrix(file4_NOS);
Data_0001PAL_NOS0920 = readmatrix(file0920_1);
Data_0001PAL_NOS_SR = readmatrix(file0920_2);
% Converged Data without scattering and H2O/CO2 absorption
file1_NOSNOA = fullfile(dir1,dir2,'OUTPUT_PLOT1PAL_NOSNOA.dat');
file2_NOSNOA = fullfile(dir1,dir2,'OUTPUT_PLOT01PAL_NOSNOA.dat');
file3_NOSNOA = fullfile(dir1,dir2,'OUTPUT_PLOT001PAL_NOSNOA.dat');
file4_NOSNOA = fullfile(dir1,dir2,'OUTPUT_PLOT0001PAL_NOSNOA.dat');
Data_1PAL_1D_NOSNOA = readmatrix(file1_NOSNOA);
Data_01PAL_1D_NOSNOA = readmatrix(file2_NOSNOA);
Data_001PAL_1D_NOSNOA = readmatrix(file3_NOSNOA);
Data_0001PAL_1D_NOSNOA = readmatrix(file4_NOSNOA);

% One step from the standard without scattering 
file1_NOS_ONESTEP = fullfile(dir1,dir3,'OUTPUT_PLOT1PAL_NOS.dat');
file2_NOS_ONESTEP = fullfile(dir1,dir3,'OUTPUT_PLOT01PAL_NOS.dat');
file3_NOS_ONESTEP = fullfile(dir1,dir3,'OUTPUT_PLOT001PAL_NOS.dat');
file4_NOS_ONESTEP = fullfile(dir1,dir3,'OUTPUT_PLOT0001PAL_NOS.dat');
Data_1PAL_1D_NOS_ONESTEP = readmatrix(file1_NOS_ONESTEP);
Data_01PAL_1D_NOS_ONESTEP = readmatrix(file2_NOS_ONESTEP);
Data_001PAL_1D_NOS_ONESTEP = readmatrix(file3_NOS_ONESTEP);
Data_0001PAL_1D_NOS_ONESTEP = readmatrix(file4_NOS_ONESTEP);
% One step from the standard without scattering and H2O/CO2 absorption
file1_NOSNOA_ONESTEP = fullfile(dir1,dir3,'OUTPUT_PLOT1PAL_NOSNOA.dat');
file2_NOSNOA_ONESTEP = fullfile(dir1,dir3,'OUTPUT_PLOT01PAL_NOSNOA.dat');
file3_NOSNOA_ONESTEP = fullfile(dir1,dir3,'OUTPUT_PLOT001PAL_NOSNOA.dat');
file4_NOSNOA_ONESTEP = fullfile(dir1,dir3,'OUTPUT_PLOT0001PAL_NOSNOA.dat');
Data_1PAL_1D_NOSNOA_ONESTEP = readmatrix(file1_NOSNOA_ONESTEP);
Data_01PAL_1D_NOSNOA_ONESTEP = readmatrix(file2_NOSNOA_ONESTEP);
Data_001PAL_1D_NOSNOA_ONESTEP = readmatrix(file3_NOSNOA_ONESTEP);
Data_0001PAL_1D_NOSNOA_ONESTEP = readmatrix(file4_NOSNOA_ONESTEP);

% read altitude profile at different O2 levels
alt = DATA_1PAL_1D(:,1)./1e5; %km

% read pressure profile at different O2 levels
Press_1PAL = DATA_1PAL_1D(:,2)./1e3; %pressure (hPa or mbar)
Press_01PAL = DATA_01PAL_1D(:,2)./1e3;
Press_001PAL = DATA_001PAL_1D(:,2)./1e3;
Press_0001PAL = DATA_0001PAL_1D(:,2)./1e3;
%% read standard converged run
% read total density profile at different O2 levels (cm^-3)
DEN_1PAL = DATA_1PAL_1D(:,3);
DEN_01PAL = DATA_01PAL_1D(:,3);
DEN_001PAL = DATA_001PAL_1D(:,3);
DEN_0001PAL = DATA_0001PAL_1D(:,3);

% O(1D) + H2O --> 2OH
R1OH_1PAL = DATA_1PAL_1D(:,4);
R1OH_01PAL = DATA_01PAL_1D(:,4);
R1OH_001PAL = DATA_001PAL_1D(:,4);
R1OH_0001PAL = DATA_0001PAL_1D(:,4);
% H2O + HV --> H + OH
R25OH_1PAL = DATA_1PAL_1D(:,5);
R25OH_01PAL = DATA_01PAL_1D(:,5);
R25OH_001PAL = DATA_001PAL_1D(:,5);
R25OH_0001PAL = DATA_0001PAL_1D(:,5);

% photolysis rate of O2 in the unit of s^-1
PO2_1PAL = DATA_1PAL_1D(:,8);
PO2_01PAL = DATA_01PAL_1D(:,8);
PO2_001PAL = DATA_001PAL_1D(:,8);
PO2_0001PAL = DATA_0001PAL_1D(:,8);
% photolysis rate of H2O in the unit of s^-1
PH2O_1PAL = DATA_1PAL_1D(:,9);
PH2O_01PAL = DATA_01PAL_1D(:,9);
PH2O_001PAL = DATA_001PAL_1D(:,9);
PH2O_0001PAL = DATA_0001PAL_1D(:,9);
% PO3D in the unit of s^-1
PO3D_1PAL = DATA_1PAL_1D(:,10);
PO3D_01PAL = DATA_01PAL_1D(:,10);
PO3D_001PAL = DATA_001PAL_1D(:,10);
PO3D_0001PAL = DATA_0001PAL_1D(:,10);
% photolysis rate of H2O2 in the unit of s^-1
PH2O2_1PAL = DATA_1PAL_1D(:,11);
PH2O2_01PAL = DATA_01PAL_1D(:,11);
PH2O2_001PAL = DATA_001PAL_1D(:,11);
PH2O2_0001PAL = DATA_0001PAL_1D(:,11);

% read mixing ratios of O3 at different O2 levels
mrO3_1PAL = DATA_1PAL_1D(:,12);
mrO3_01PAL = DATA_01PAL_1D(:,12);
mrO3_001PAL = DATA_001PAL_1D(:,12);
mrO3_0001PAL = DATA_0001PAL_1D(:,12);

% read mixing ratios of H2O at different O2 levels
mrH2O_1PAL = DATA_1PAL_1D(:,13);
mrH2O_01PAL = DATA_01PAL_1D(:,13);
mrH2O_001PAL = DATA_001PAL_1D(:,13);
mrH2O_0001PAL = DATA_0001PAL_1D(:,13);
% read number density of O2 at different O2 levels
NumO2_1PAL = DATA_1PAL_1D(:,14);
NumO2_01PAL = DATA_01PAL_1D(:,14);
NumO2_001PAL = DATA_001PAL_1D(:,14);
NumO2_0001PAL = DATA_0001PAL_1D(:,14);
% read number density of OH at different O2 levels
NumOH_1PAL = DATA_1PAL_1D(:,15);
NumOH_01PAL = DATA_01PAL_1D(:,15);
NumOH_001PAL = DATA_001PAL_1D(:,15);
NumOH_0001PAL = DATA_0001PAL_1D(:,15);
% read number density of O1D at different O2 levels
NumO1D_1PAL = DATA_1PAL_1D(:,16);
NumO1D_01PAL = DATA_01PAL_1D(:,16);
NumO1D_001PAL = DATA_001PAL_1D(:,16);
NumO1D_0001PAL = DATA_0001PAL_1D(:,16);
% calculate number density of O3: DENSITY*MXIING RATIOS for 1D
NumO3_1PAL = DEN_1PAL.*mrO3_1PAL;
NumO3_01PAL = DEN_01PAL.*mrO3_01PAL;
NumO3_001PAL = DEN_001PAL.*mrO3_001PAL;
NumO3_0001PAL = DEN_0001PAL.*mrO3_0001PAL;

%% Effect of scattering and absorbers on the photolysis rates of O2 and H2O in the 1D 
% read ONESTEP DATA
% only without scattering in 1D
% 
PO2_NOS_1PAL = Data_1PAL_1D_NOS_ONESTEP(:,8);
PO2_NOS_01PAL = Data_01PAL_1D_NOS_ONESTEP(:,8);
PO2_NOS_001PAL = Data_001PAL_1D_NOS_ONESTEP(:,8);
PO2_NOS_0001PAL = Data_0001PAL_1D_NOS_ONESTEP(:,8);
% PO3D in the unit of s^-1
PO3D_NOS_1PAL = Data_1PAL_1D_NOS_ONESTEP(:,10);
PO3D_NOS_01PAL = Data_01PAL_1D_NOS_ONESTEP(:,10);
PO3D_NOS_001PAL = Data_001PAL_1D_NOS_ONESTEP(:,10);
PO3D_NOS_0001PAL = Data_0001PAL_1D_NOS_ONESTEP(:,10);
% photolysis rate of H2O in the unit of s^-1
PH2O_NOS_1PAL = Data_1PAL_1D_NOS_ONESTEP(:,9);
PH2O_NOS_01PAL = Data_01PAL_1D_NOS_ONESTEP(:,9);
PH2O_NOS_001PAL = Data_001PAL_1D_NOS_ONESTEP(:,9);
PH2O_NOS_0001PAL = Data_0001PAL_1D_NOS_ONESTEP(:,9);
% photolysis rate of H2O2 in the unit of s^-1
PH2O2_NOS_1PAL = Data_1PAL_1D_NOS_ONESTEP(:,11);
PH2O2_NOS_01PAL = Data_01PAL_1D_NOS_ONESTEP(:,11);
PH2O2_NOS_001PAL = Data_001PAL_1D_NOS_ONESTEP(:,11);
PH2O2_NOS_0001PAL = Data_0001PAL_1D_NOS_ONESTEP(:,11);

% without scattering/absorbers in 1D
PO2_NO_1PAL = Data_1PAL_1D_NOSNOA_ONESTEP(:,8);
PO2_NO_01PAL = Data_01PAL_1D_NOSNOA_ONESTEP(:,8);
PO2_NO_001PAL = Data_001PAL_1D_NOSNOA_ONESTEP(:,8);
PO2_NO_0001PAL = Data_0001PAL_1D_NOSNOA_ONESTEP(:,8);

% photolysis rate of H2O in the unit of s^-1
PH2O_NO_1PAL = Data_1PAL_1D_NOSNOA_ONESTEP(:,9);
PH2O_NO_01PAL = Data_01PAL_1D_NOSNOA_ONESTEP(:,9);
PH2O_NO_001PAL = Data_001PAL_1D_NOSNOA_ONESTEP(:,9);
PH2O_NO_0001PAL = Data_0001PAL_1D_NOSNOA_ONESTEP(:,9);

% PO3D in the unit of s^-1
PO3D_NO_1PAL = Data_1PAL_1D_NOSNOA_ONESTEP(:,10);
PO3D_NO_01PAL = Data_01PAL_1D_NOSNOA_ONESTEP(:,10);
PO3D_NO_001PAL = Data_001PAL_1D_NOSNOA_ONESTEP(:,10);
PO3D_NO_0001PAL = Data_0001PAL_1D_NOSNOA_ONESTEP(:,10);

% photolysis rate of H2O2 in the unit of s^-1
PH2O2_NO_1PAL = Data_1PAL_1D_NOSNOA_ONESTEP(:,11);
PH2O2_NO_01PAL = Data_01PAL_1D_NOSNOA_ONESTEP(:,11);
PH2O2_NO_001PAL = Data_001PAL_1D_NOSNOA_ONESTEP(:,11);
PH2O2_NO_0001PAL = Data_0001PAL_1D_NOSNOA_ONESTEP(:,11);
%% Number density without scattering/absorbers
% read Converged DATA
% only without scattering
R1OH_NOS_1PAL = Data_1PAL_1D_NOS(:,4);
R1OH_NOS_01PAL = Data_01PAL_1D_NOS(:,4);
R1OH_NOS_001PAL = Data_001PAL_1D_NOS(:,4);
R1OH_NOS_0001PAL = Data_0001PAL_1D_NOS(:,4);

R25OH_NOS_1PAL = Data_1PAL_1D_NOS(:,5);
R25OH_NOS_01PAL = Data_01PAL_1D_NOS(:,5);
R25OH_NOS_001PAL = Data_001PAL_1D_NOS(:,5);
R25OH_NOS_0001PAL = Data_0001PAL_1D_NOS(:,5);
R25OH_NOS_0001PAL0920 = Data_0001PAL_NOS0920(:,5);
R25OH_NOS_0001PAL_SR = Data_0001PAL_NOS_SR(:,5);
PO2_NOS_0001PAL0920 = Data_0001PAL_NOS0920(:,8);
NumO2_NOS_0001PAL0920 = Data_0001PAL_NOS0920(:,14);
PO2_NOS_0001PAL_SR = Data_0001PAL_NOS_SR(:,8);
NumO2_NOS_0001PAL_SR = Data_0001PAL_NOS_SR(:,14);

NumO2_NOS_1PAL = Data_1PAL_1D_NOS(:,14);
NumO2_NOS_01PAL = Data_01PAL_1D_NOS(:,14);
NumO2_NOS_001PAL = Data_001PAL_1D_NOS(:,14);
NumO2_NOS_0001PAL = Data_0001PAL_1D_NOS(:,14);
NumOH_NOS_1PAL = Data_1PAL_1D_NOS(:,15);
NumOH_NOS_01PAL = Data_01PAL_1D_NOS(:,15);
NumOH_NOS_001PAL = Data_001PAL_1D_NOS(:,15);
NumOH_NOS_0001PAL = Data_0001PAL_1D_NOS(:,15);
% without scattering and absorbers
R1OH_NO_1PAL = Data_1PAL_1D_NOSNOA(:,4);
R1OH_NO_01PAL = Data_01PAL_1D_NOSNOA(:,4);
R1OH_NO_001PAL = Data_001PAL_1D_NOSNOA(:,4);
R1OH_NO_0001PAL = Data_0001PAL_1D_NOSNOA(:,4);

R25OH_NO_1PAL = Data_1PAL_1D_NOSNOA(:,5);
R25OH_NO_01PAL = Data_01PAL_1D_NOSNOA(:,5);
R25OH_NO_001PAL = Data_001PAL_1D_NOSNOA(:,5);
R25OH_NO_0001PAL = Data_0001PAL_1D_NOSNOA(:,5);

NumO2_NO_1PAL = Data_1PAL_1D_NOSNOA(:,14);
NumO2_NO_01PAL = Data_01PAL_1D_NOSNOA(:,14);
NumO2_NO_001PAL = Data_001PAL_1D_NOSNOA(:,14);
NumO2_NO_0001PAL = Data_0001PAL_1D_NOSNOA(:,14);
NumOH_NO_1PAL = Data_1PAL_1D_NOSNOA(:,15);
NumOH_NO_01PAL = Data_01PAL_1D_NOSNOA(:,15);
NumOH_NO_001PAL = Data_001PAL_1D_NOSNOA(:,15);
NumOH_NO_0001PAL = Data_0001PAL_1D_NOSNOA(:,15);
%% ********** Production rate of Ox in 1d **********
% standar run
ProOx_1PAL_1D = PO2_1PAL.*NumO2_1PAL;
ProOx_01PAL_1D = PO2_01PAL.*NumO2_01PAL;
ProOx_001PAL_1D = PO2_001PAL.*NumO2_001PAL;
ProOx_0001PAL_1D = PO2_0001PAL.*NumO2_0001PAL;
% without scattering/absorbers in 1d
ProOx_NO_1PAL = PO2_NO_1PAL.*NumO2_NO_1PAL;
ProOx_NO_01PAL = PO2_NO_01PAL.*NumO2_NO_01PAL;
ProOx_NO_001PAL = PO2_NO_001PAL.*NumO2_NO_001PAL;
ProOx_NO_0001PAL = PO2_NO_0001PAL.*NumO2_NO_0001PAL;
% ONLY without scattering in 1d
ProOx_NOS_1PAL = PO2_NOS_1PAL.*NumO2_NOS_1PAL;
ProOx_NOS_01PAL = PO2_NOS_01PAL.*NumO2_NOS_01PAL;
ProOx_NOS_001PAL = PO2_NOS_001PAL.*NumO2_NOS_001PAL;
ProOx_NOS_0001PAL = PO2_NOS_0001PAL.*NumO2_NOS_0001PAL;
ProOx_NOS_0001PAL0920 = PO2_NOS_0001PAL0920.*NumO2_NOS_0001PAL0920;
ProOx_NOS_0001PAL_SR = PO2_NOS_0001PAL_SR.*NumO2_NOS_0001PAL_SR;

%% **************** read data from WACCM3D ****************
DATA_NO3_1PAL_WACCM3D = readmatrix('PI_j_mr_T_Z3_data_IC.csv');
DATA_NO3_01PAL_WACCM3D = readmatrix('Ten_j_mr_T_Z3_data_IC.csv');
DATA_NO3_001PAL_WACCM3D = readmatrix('One_j_mr_T_Z3_data_IC.csv');
DATA_NO3_0001PAL_WACCM3D = readmatrix('Zero1_j_mr_T_Z3_data_IC.csv');
% ********** read each column in WACCM3D ************
% pressure (hPa)
TotalP_1PAL = flip(DATA_NO3_1PAL_WACCM3D(:,3));
TotalP_01PAL = flip(DATA_NO3_01PAL_WACCM3D(:,3));
TotalP_001PAL = flip(DATA_NO3_001PAL_WACCM3D(:,3));
TotalP_0001PAL = flip(DATA_NO3_0001PAL_WACCM3D(:,3));

% mixing ratio of O1D
mrO1D_1PAL_WACCM3D = flip(DATA_NO3_1PAL_WACCM3D(:,9));
mrO1D_01PAL_WACCM3D = flip(DATA_NO3_01PAL_WACCM3D(:,9));
mrO1D_001PAL_WACCM3D = flip(DATA_NO3_001PAL_WACCM3D(:,9));
mrO1D_0001PAL_WACCM3D = flip(DATA_NO3_0001PAL_WACCM3D(:,9));
% mixing ratio of O2
mrO2_1PAL_WACCM3D = flip(DATA_NO3_1PAL_WACCM3D(:,11));
mrO2_01PAL_WACCM3D = flip(DATA_NO3_01PAL_WACCM3D(:,11));
mrO2_001PAL_WACCM3D = flip(DATA_NO3_001PAL_WACCM3D(:,11));
mrO2_0001PAL_WACCM3D = flip(DATA_NO3_0001PAL_WACCM3D(:,11));
% mixing ratio of O3
mrO3_1PAL_WACCM3D = flip(DATA_NO3_1PAL_WACCM3D(:,12));
mrO3_01PAL_WACCM3D = flip(DATA_NO3_01PAL_WACCM3D(:,12));
mrO3_001PAL_WACCM3D = flip(DATA_NO3_001PAL_WACCM3D(:,12));
mrO3_0001PAL_WACCM3D = flip(DATA_NO3_0001PAL_WACCM3D(:,12));
% mixing ratios of H2O
mrH2O_1PAL_WACCM3D = flip(DATA_NO3_1PAL_WACCM3D(:,13));
mrH2O_01PAL_WACCM3D = flip(DATA_NO3_01PAL_WACCM3D(:,13));
mrH2O_001PAL_WACCM3D = flip(DATA_NO3_001PAL_WACCM3D(:,13));
mrH2O_0001PAL_WACCM3D = flip(DATA_NO3_0001PAL_WACCM3D(:,13));

% photolysis rate of O2 (s^-1)
PO2_1PAL_WACCM3D = flip(DATA_NO3_1PAL_WACCM3D(:,15));
PO2_01PAL_WACCM3D = flip(DATA_NO3_01PAL_WACCM3D(:,15));
PO2_001PAL_WACCM3D = flip(DATA_NO3_001PAL_WACCM3D(:,15));
PO2_0001PAL_WACCM3D = flip(DATA_NO3_0001PAL_WACCM3D(:,15));
% photolysis rate of H2O (s^-1)
PH2O_1PAL_WACCM3D = flip(DATA_NO3_1PAL_WACCM3D(:,17));
PH2O_01PAL_WACCM3D = flip(DATA_NO3_01PAL_WACCM3D(:,17));
PH2O_001PAL_WACCM3D = flip(DATA_NO3_001PAL_WACCM3D(:,17));
PH2O_0001PAL_WACCM3D = flip(DATA_NO3_0001PAL_WACCM3D(:,17));
%% mixing ratios of CH4 (may not be used)
fCH4_1PAL_WACCM3D = flip(DATA_NO3_1PAL_WACCM3D(:,14));
fCH4_01PAL_WACCM3D = flip(DATA_NO3_01PAL_WACCM3D(:,14));
fCH4_001PAL_WACCM3D = flip(DATA_NO3_001PAL_WACCM3D(:,14));
fCH4_0001PAL_WACCM3D = flip(DATA_NO3_0001PAL_WACCM3D(:,14));
%% calculate density in WACCM3D with its temperature and pressure
kB = 1.38e-23; % J K^-1
T_1PAL_WACCM3D = flip(DATA_NO3_1PAL_WACCM3D(:,2));
T_01PAL_WACCM3D = flip(DATA_NO3_01PAL_WACCM3D(:,2));
T_001PAL_WACCM3D = flip(DATA_NO3_001PAL_WACCM3D(:,2));
T_0001PAL_WACCM3D = flip(DATA_NO3_0001PAL_WACCM3D(:,2));
DEN_1PAL_WACCM3D = TotalP_1PAL./(kB*T_1PAL_WACCM3D)*100*1e-6;
DEN_01PAL_WACCM3D = TotalP_01PAL./(kB*T_01PAL_WACCM3D)*100*1e-6;
DEN_001PAL_WACCM3D = TotalP_001PAL./(kB*T_001PAL_WACCM3D)*100*1e-6;
DEN_0001PAL_WACCM3D = TotalP_0001PAL./(kB*T_0001PAL_WACCM3D)*100*1e-6;

% calculate the number density (cm^-3) in WACCM3D
NumH2O_1PAL_WACCM3D = DEN_1PAL_WACCM3D.*mrH2O_1PAL_WACCM3D;
NumH2O_01PAL_WACCM3D = DEN_01PAL_WACCM3D.*mrH2O_01PAL_WACCM3D;
NumH2O_001PAL_WACCM3D = DEN_001PAL_WACCM3D.*mrH2O_001PAL_WACCM3D;
NumH2O_0001PAL_WACCM3D = DEN_0001PAL_WACCM3D.*mrH2O_0001PAL_WACCM3D;
% calculate number density of O3 in WACCM6
NumO3_1PAL_WACCM3D = DEN_1PAL_WACCM3D.*mrO3_1PAL_WACCM3D;
NumO3_01PAL_WACCM3D = DEN_01PAL_WACCM3D.*mrO3_01PAL_WACCM3D;
NumO3_001PAL_WACCM3D = DEN_001PAL_WACCM3D.*mrO3_001PAL_WACCM3D;
NumO3_0001PAL_WACCM3D = DEN_0001PAL_WACCM3D.*mrO3_0001PAL_WACCM3D;
% Number density of O2 in WACCM
NO2_1PAL_WACCM3D = mrO2_1PAL_WACCM3D.*DEN_1PAL_WACCM3D;
NO2_01PAL_WACCM3D = mrO2_01PAL_WACCM3D.*DEN_01PAL_WACCM3D;
NO2_001PAL_WACCM3D = mrO2_001PAL_WACCM3D.*DEN_001PAL_WACCM3D;
NO2_0001PAL_WACCM3D = mrO2_0001PAL_WACCM3D.*DEN_0001PAL_WACCM3D;
%%
NumO1D_1PAL_WACCM3D = DEN_1PAL_WACCM3D.*mrO1D_1PAL_WACCM3D;
NumO1D_01PAL_WACCM3D = DEN_01PAL_WACCM3D.*mrO1D_01PAL_WACCM3D;
NumO1D_001PAL_WACCM3D = DEN_001PAL_WACCM3D.*mrO1D_001PAL_WACCM3D;
NumO1D_0001PAL_WACCM3D = DEN_0001PAL_WACCM3D.*mrO1D_0001PAL_WACCM3D;
%% ********** Production rate of Ox in WACCM3D **********
ProOx_1PAL_WACCM3D = PO2_1PAL_WACCM3D.*NO2_1PAL_WACCM3D;
ProOx_01PAL_WACCM3D = PO2_01PAL_WACCM3D.*NO2_01PAL_WACCM3D;
ProOx_001PAL_WACCM3D = PO2_001PAL_WACCM3D.*NO2_001PAL_WACCM3D;
ProOx_0001PAL_WACCM3D = PO2_0001PAL_WACCM3D.*NO2_0001PAL_WACCM3D;
%% OH number density in the WACCM6 3D
DATA_1PAL_WACCM3D = readmatrix('PI_OH_density.csv');
DATA_01PAL_WACCM3D = readmatrix('Ten_OH_density.csv');
DATA_001PAL_WACCM3D = readmatrix('One_OH_density.csv');
DATA_0001PAL_WACCM3D = readmatrix('Zero1_OH_density.csv');
NOH_1PAL_WACCM3D = DATA_1PAL_WACCM3D(:,4);
NOH_01PAL_WACCM3D = DATA_01PAL_WACCM3D(:,4);
NOH_001PAL_WACCM3D = DATA_001PAL_WACCM3D(:,4);
NOH_0001PAL_WACCM3D = DATA_0001PAL_WACCM3D(:,4);
NP_1PAL_WACCM3D = DATA_1PAL_WACCM3D(:,2);
NP_01PAL_WACCM3D = DATA_01PAL_WACCM3D(:,2);
NP_001PAL_WACCM3D = DATA_001PAL_WACCM3D(:,2);
NP_0001PAL_WACCM3D = DATA_0001PAL_WACCM3D(:,2);
%% RATE OF OH SOURCE FROM H2O PHOTOLYSIS IN WACCM 3D
% OH from H2O photolysis
ROH1PAL_WACCM3D = NumH2O_1PAL_WACCM3D.*PH2O_1PAL_WACCM3D;
ROH01PAL_WACCM3D = NumH2O_01PAL_WACCM3D.*PH2O_01PAL_WACCM3D;
ROH001PAL_WACCM3D = NumH2O_001PAL_WACCM3D.*PH2O_001PAL_WACCM3D;
ROH0001PAL_WACCM3D = NumH2O_0001PAL_WACCM3D.*PH2O_0001PAL_WACCM3D;
% OH from O1D + H2O
RateOH1PAL_WACCM3D = readmatrix('PI_O1D_H2O.csv');
RateOH01PAL_WACCM3D = readmatrix('Ten_O1D_H2O.csv');
RateOH001PAL_WACCM3D = readmatrix('One_O1D_H2O.csv');
RateOH0001PAL_WACCM3D = readmatrix('Zero1_O1D_H2O.csv');
% calculate density in WACCM3D with its temperature and pressure
TotalP_1PAL_3D = flip(RateOH1PAL_WACCM3D(:,3));
TotalP_01PAL_3D = flip(RateOH01PAL_WACCM3D(:,3));
TotalP_001PAL_3D = flip(RateOH001PAL_WACCM3D(:,3));
TotalP_0001PAL_3D = flip(RateOH0001PAL_WACCM3D(:,3));
kB = 1.38e-23; % J K^-1
Temp_1PAL_3D = flip(RateOH1PAL_WACCM3D(:,2));
Temp_01PAL_3D = flip(RateOH01PAL_WACCM3D(:,2));
Temp_001PAL_3D = flip(RateOH001PAL_WACCM3D(:,2));
Temp_0001PAL_3D = flip(RateOH0001PAL_WACCM3D(:,2));

DEN_1PAL_3D = TotalP_1PAL_3D./(kB*Temp_1PAL_3D)*100*1e-6;
DEN_01PAL_3D = TotalP_01PAL_3D./(kB*Temp_01PAL_3D)*100*1e-6;
DEN_001PAL_3D = TotalP_001PAL_3D./(kB*Temp_001PAL_3D)*100*1e-6;
DEN_0001PAL_3D = TotalP_0001PAL_3D./(kB*Temp_0001PAL_3D)*100*1e-6;

mrH2O_1PAL_3D = flip(RateOH1PAL_WACCM3D(:,5));
mrH2O_01PAL_3D = flip(RateOH01PAL_WACCM3D(:,5));
mrH2O_001PAL_3D = flip(RateOH001PAL_WACCM3D(:,5));
mrH2O_0001PAL_3D = flip(RateOH0001PAL_WACCM3D(:,5));

mrO1D_1PAL_3D = flip(RateOH1PAL_WACCM3D(:,6));
mrO1D_01PAL_3D = flip(RateOH01PAL_WACCM3D(:,6));
mrO1D_001PAL_3D = flip(RateOH001PAL_WACCM3D(:,6));
mrO1D_0001PAL_3D = flip(RateOH0001PAL_WACCM3D(:,6));

ROH_1PAL_3D = flip(RateOH1PAL_WACCM3D(:,7));
ROH_01PAL_3D = flip(RateOH01PAL_WACCM3D(:,7));
ROH_001PAL_3D = flip(RateOH001PAL_WACCM3D(:,7));
ROH_0001PAL_3D = flip(RateOH0001PAL_WACCM3D(:,7));

%Number density of H2O and O1D
NumH2O_1PAL_3D = mrH2O_1PAL_3D.*DEN_1PAL_3D;
NumH2O_01PAL_3D = mrH2O_01PAL_3D.*DEN_01PAL_3D;
NumH2O_001PAL_3D = mrH2O_001PAL_3D.*DEN_001PAL_3D;
NumH2O_0001PAL_3D = mrH2O_0001PAL_3D.*DEN_0001PAL_3D;

NumO1D_1PAL_3D = mrO1D_1PAL_3D.*DEN_1PAL_3D;
NumO1D_01PAL_3D = mrO1D_01PAL_3D.*DEN_01PAL_3D;
NumO1D_001PAL_3D = mrO1D_001PAL_3D.*DEN_001PAL_3D;
NumO1D_0001PAL_3D = mrO1D_0001PAL_3D.*DEN_0001PAL_3D;

% OH production rate
ProOH_1PAL_3D = ROH_1PAL_3D.*NumH2O_1PAL_3D.*NumO1D_1PAL_3D;
ProOH_01PAL_3D = ROH_01PAL_3D.*NumH2O_01PAL_3D.*NumO1D_01PAL_3D;
ProOH_001PAL_3D = ROH_001PAL_3D.*NumH2O_001PAL_3D.*NumO1D_001PAL_3D;
ProOH_0001PAL_3D = ROH_0001PAL_3D.*NumH2O_0001PAL_3D.*NumO1D_0001PAL_3D;
%% New simulations from WACCM6 with H2O/CO2 absorption
DATA_01PAL_WACCM1206 = readmatrix('Updated_run_data_0.1_PAL_IC_1206.csv');
%DATA_001PAL_WACCM1206 = readmatrix('Updated_run_data_0.01_PAL_IC_1206.csv');
DATA_001PAL_WACCM1206 = readmatrix('Updated_run_data_0.01_PAL_SRB_abs_and_LBC_IC1209.csv');
DATA_0001PAL_WACCM1206 = readmatrix('Updated_run_data_0.001_PAL_IC_1206.csv');
Press_01PAL_WACCM1206 = flip(DATA_01PAL_WACCM1206(:,3));
Press_001PAL_WACCM1206 = flip(DATA_001PAL_WACCM1206(:,3));
Press_0001PAL_WACCM1206 = flip(DATA_0001PAL_WACCM1206(:,3));
mrH2O_01PAL_WACCM1206 = flip(DATA_01PAL_WACCM1206(:,8));
mrH2O_001PAL_WACCM1206 = flip(DATA_001PAL_WACCM1206(:,8));
mrH2O_0001PAL_WACCM1206 = flip(DATA_0001PAL_WACCM1206(:,8));
mrOH_01PAL_WACCM1206 = flip(DATA_01PAL_WACCM1206(:,13));
mrOH_001PAL_WACCM1206 = flip(DATA_001PAL_WACCM1206(:,13));
mrOH_0001PAL_WACCM1206 = flip(DATA_0001PAL_WACCM1206(:,13));
TotalNum_01PAL_WACCM1206 = flip(DATA_01PAL_WACCM1206(:,15));
TotalNum_001PAL_WACCM1206 = flip(DATA_001PAL_WACCM1206(:,15));
TotalNum_0001PAL_WACCM1206 = flip(DATA_0001PAL_WACCM1206(:,15));
PH2O_01PAL_WACCM1206 = flip(DATA_01PAL_WACCM1206(:,10));
PH2O_001PAL_WACCM1206 = flip(DATA_001PAL_WACCM1206(:,10));
PH2O_0001PAL_WACCM1206 = flip(DATA_0001PAL_WACCM1206(:,10));

NumH2O_01PAL_WACCM1206 = mrH2O_01PAL_WACCM1206.*TotalNum_01PAL_WACCM1206;
NumH2O_001PAL_WACCM1206 = mrH2O_001PAL_WACCM1206.*TotalNum_001PAL_WACCM1206;
NumH2O_0001PAL_WACCM1206 = mrH2O_0001PAL_WACCM1206.*TotalNum_0001PAL_WACCM1206;
NumOH_01PAL_WACCM1206 = mrOH_01PAL_WACCM1206.*TotalNum_01PAL_WACCM1206;
NumOH_001PAL_WACCM1206 = mrOH_001PAL_WACCM1206.*TotalNum_001PAL_WACCM1206;
NumOH_0001PAL_WACCM1206 = mrOH_0001PAL_WACCM1206.*TotalNum_0001PAL_WACCM1206;
RateOH_01PAL_WACCM1206 = NumH2O_01PAL_WACCM1206.*PH2O_01PAL_WACCM1206;
RateOH_001PAL_WACCM1206 = NumH2O_001PAL_WACCM1206.*PH2O_001PAL_WACCM1206;
RateOH_0001PAL_WACCM1206 = NumH2O_0001PAL_WACCM1206.*PH2O_0001PAL_WACCM1206;

%% **************** read data from ROCKE3D ****************
%ncdisp('profiles.nc') check information in the data
Press_ROCK3D = ncread('profiles.nc','pressure'); 
pO2_ROCK3D = ncread('profiles.nc','pO2');
Temp_ROCK3D = ncread('profiles.nc','TempL');
Ox_ROCK3D = ncread('profiles.nc','Ox');
NumOH_ROCK3D = ncread('profiles.nc','OH_conc');
H2O_ROCK3D = ncread('profiles.nc','Water');
JO1D_ROCK3D = ncread('profiles.nc','JO1D');
fCH4_ROCK3D = ncread('profiles.nc','CH4');
% ***** calculate denisty in ROCKE3D with its temperature and pressure ****
DEN_1PAL_ROCKE3D = Press_ROCK3D./(kB*Temp_ROCK3D(:,1))*100*1e-6;
DEN_01PAL_ROCKE3D = Press_ROCK3D./(kB*Temp_ROCK3D(:,2))*100*1e-6;
DEN_001PAL_ROCKE3D = Press_ROCK3D./(kB*Temp_ROCK3D(:,3))*100*1e-6;
DEN_0001PAL_ROCKE3D = Press_ROCK3D./(kB*Temp_ROCK3D(:,4))*100*1e-6;
% calculate number density of Ox in ROCKE3D
NumOx_1PAL_ROCKE3D = DEN_1PAL_ROCKE3D.*Ox_ROCK3D(:,1);
NumOx_01PAL_ROCKE3D = DEN_01PAL_ROCKE3D.*Ox_ROCK3D(:,2);
NumOx_001PAL_ROCKE3D = DEN_001PAL_ROCKE3D.*Ox_ROCK3D(:,3);
NumOx_0001PAL_ROCKE3D = DEN_0001PAL_ROCKE3D.*Ox_ROCK3D(:,4);

%% ***************** read data from Kinetics 1D ****************
%ncdisp('CaltechJPL_prodOx_4aoshuang.nc')
Press_Kinetics = ncread('CaltechJPL_prodOx_4aoshuang.nc','p');
ProOx_1PAL_Kinetics = ncread('CaltechJPL_prodOx_4aoshuang.nc','pox_1pal');
ProOx_01PAL_Kinetics = ncread('CaltechJPL_prodOx_4aoshuang.nc','pox_0p1pal');
ProOx_001PAL_Kinetics = ncread('CaltechJPL_prodOx_4aoshuang.nc','pox_0p01pal');
ProOx_0001PAL_Kinetics = ncread('CaltechJPL_prodOx_4aoshuang.nc','pox_0p001pal');

%% Figure 2. Temperature and eddy diffusion plot profile
% ********* TEMPERATURE PROFILES FOR DIFFERENT O2 LEVELS IN 1D ***********
Temp1D = readmatrix('Temp_1DPHOTO.csv'); 
Temp = Temp1D(:,2); % 1PAL pO2
Temp1 = Temp1D(:,3); % 0.1PAL pO2
Temp2 = Temp1D(:,4); % <0.1PAL pO2
% ********* TEMPERATURE PROFILES FROM WACCM6 3D ********************
DATA_Temp_1PAL_WACCM3D = readmatrix('PI_jO2_T_SZA.csv');
DATA_Temp_01PAL_WACCM3D = readmatrix('Ten_jO2_T_SZA.csv');
DATA_Temp_001PAL_WACCM3D = readmatrix('One_jO2_T_SZA.csv');
% Pressure in WACCM6
Press_1PAL_WACCM3D = flip(DATA_Temp_1PAL_WACCM3D(:,2));
Press_01PAL_WACCM3D = flip(DATA_Temp_01PAL_WACCM3D(:,2));
Press_001PAL_WACCM3D = flip(DATA_Temp_001PAL_WACCM3D(:,2));
% Global mean temperature
Temp_1PAL_WACCM3D = flip(DATA_Temp_1PAL_WACCM3D(:,5));
Temp_01PAL_WACCM3D = flip(DATA_Temp_01PAL_WACCM3D(:,5));
Temp_001PAL_WACCM3D = flip(DATA_Temp_001PAL_WACCM3D(:,5));

% ********* EDDY DIFFUSION IN 1D ***********
eddy = readmatrix('eddy_1D');

% Plot temperature profiles against pressure
figure
subplot(121)
plot(Temp,Press_1PAL,'Color',[0.4940 0.1840 0.5560],'LineWidth',2)
hold on
plot(Temp1,Press_01PAL,'Color',[0.85,0.33,0.10],'LineWidth',2)
plot(Temp2,Press_001PAL,'Color',[0.93,0.69,0.13],'LineWidth',2)
plot(Temp_1PAL_WACCM3D,Press_1PAL_WACCM3D,'-.','Color',[0.4940 0.1840 0.5560],'LineWidth',2)
plot(Temp_01PAL_WACCM3D,Press_01PAL_WACCM3D,'-.','Color',[0.85,0.33,0.10],'LineWidth',2)
plot(Temp_001PAL_WACCM3D,Press_001PAL_WACCM3D,'-.','Color',[0.93,0.69,0.13],'LineWidth',2)
ylim([1e-3 1e3]);
xlabel('Temperature (K)','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
legend('1PAL pO_2 in 1D','0.1PAL pO_2 in 1D','< 0.1PAL pO_2 in 1D','1 PAL pO_2 in WACCM6','0.1 PAL pO_2 in 3D','0.01 PAL pO_2 in 3D','fontsize',14);
set(gca,'LineWidth',2,'FontSize',20);
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

subplot(122)
plot(eddy(:,1),Press_1PAL,'Color',[0 0.45 0.74],'LineWidth',2)
hold on
plot(eddy(:,3),Press_1PAL,'-.','Color',[0 0.45 0.74],'LineWidth',2)
plot(eddy(:,4),Press_1PAL,'--','Color',[0 0.45 0.74],'LineWidth',2)
xlabel('Eddy diffusion (cm^{-2} s^{-1})','Fontweight','bold');
ylabel('Pressure (hPa)','Fontweight','bold');
ylim([1e-3 1e3]);
legend('Standard eddy diffusion in the 1-D','Test 1 with higher eddy diffusion in the 1-D','Test 2 with higher eddy diffusion in the 1-D','Fontsize',14);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
set(gca,'LineWidth',2,'FontSize',20);
hold off

%% Figure 3. Mixing ratio of H2O in three models
figure
subplot(221)
plot(mrH2O_1PAL,Press_1PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(mrH2O_1PAL_WACCM3D,TotalP_1PAL,'Color',[0.00,0.45,0.74],'LineWidth',2)
plot(H2O_ROCK3D(:,1),Press_ROCK3D,'Color',[0.93,0.69,0.13],'LineWidth',2)
xlabel('H_2O mixing ratio at 1 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
legend('1D','WACCM6 3D','ROCKE 3D','Fontsize',14);
xlim([1e-8 1e-2]);
ylim([1e-2 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

subplot(222)
plot(mrH2O_01PAL,Press_01PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(mrH2O_01PAL_WACCM3D,TotalP_01PAL,'Color',[0.00,0.45,0.74],'LineWidth',2)
plot(H2O_ROCK3D(:,2),Press_ROCK3D,'Color',[0.93,0.69,0.13],'LineWidth',2)
xlabel('H_2O mixing ratio at 0.1 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
xlim([1e-8 1e-2]);
ylim([1e-2 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

subplot(223)
plot(mrH2O_001PAL,Press_001PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(mrH2O_001PAL_WACCM3D,TotalP_001PAL,'Color',[0.00,0.45,0.74],'LineWidth',2)
plot(H2O_ROCK3D(:,3),Press_ROCK3D,'Color',[0.93,0.69,0.13],'LineWidth',2)
xlabel('H_2O mixing ratio at 0.01 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
xlim([1e-8 1e-2]);
ylim([1e-2 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

subplot(224)
plot(mrH2O_0001PAL,Press_0001PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(mrH2O_0001PAL_WACCM3D,TotalP_0001PAL,'Color',[0.00,0.45,0.74],'LineWidth',2)
plot(H2O_ROCK3D(:,4),Press_ROCK3D,'Color',[0.93,0.69,0.13],'LineWidth',2)
xlabel('H_2O mixing ratio at 0.0001 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
xlim([1e-8 1e-2]);
ylim([1e-2 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off
%% Figure 4:linear scale for comparing NumO3 in 1D and WACCM 3D at different O2 levels
figure
subplot(221)
plot(NumO3_1PAL./1e11,Press_1PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(NumO3_1PAL_WACCM3D./1e11,TotalP_1PAL,'Color',[0.00,0.45,0.74],'LineWidth',2)
plot(NumOx_1PAL_ROCKE3D./1e11,Press_ROCK3D,'Color',[0.93,0.69,0.13],'LineWidth',2)
xlabel('Number density of O_3 (10^{11} cm^{-3}) at 1 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
ylim([1 1e3]);
legend('1D','WACCM6 3D','ROCKE 3D');
set(gca,'YScale','log');
set(gca,'Ydir','reverse');
hold off

subplot(222)
plot(NumO3_01PAL./1e11,Press_01PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(NumO3_01PAL_WACCM3D./1e11,TotalP_01PAL,'Color',[0.00,0.45,0.74],'LineWidth',2)
plot(NumOx_01PAL_ROCKE3D./1e11,Press_ROCK3D,'Color',[0.93,0.69,0.13],'LineWidth',2)
xlabel('Number density of O_3 (10^{11} cm^{-3}) at 0.1 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
ylim([1 1e3]);
set(gca,'YScale','log');
set(gca,'Ydir','reverse');
hold off

subplot(223)
plot(NumO3_001PAL./1e11,Press_001PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(NumO3_001PAL_WACCM3D./1e11,TotalP_001PAL,'Color',[0.00,0.45,0.74],'LineWidth',2)
plot(NumOx_001PAL_ROCKE3D./1e11,Press_ROCK3D,'Color',[0.93,0.69,0.13],'LineWidth',2)
xlabel('Number density of O_3 (10^{11} cm^{-3}) at 0.01 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
ylim([1 1e3]);
set(gca,'YScale','log');
set(gca,'Ydir','reverse');
hold off

subplot(224)
plot(NumO3_0001PAL./1e11,Press_0001PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(NumO3_0001PAL_WACCM3D./1e11,TotalP_0001PAL,'Color',[0.00,0.45,0.74],'LineWidth',2)
plot(NumOx_0001PAL_ROCKE3D./1e11,Press_ROCK3D,'Color',[0.93,0.69,0.13],'LineWidth',2)
xlabel('Number density of O_3 (10^{11} cm^{-3}) at 0.001 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
ylim([1 1e3]);
set(gca,'YScale','log');
set(gca,'Ydir','reverse');
hold off
%% Figure 5. COMPARISON WITH FERNANDEZ ET AL., 2006 LINE BY LINE CALCULATION
% read PO2 data in the 1-d model and Fernandez et al.'s LBL
Data_fig5_1D = readmatrix('Figure5_1D.csv');
Z = Data_fig5_1D(:,2);
Press_1D = interp1(alt,Press_1PAL,Z./1E5);
PO2_0SZA = Data_fig5_1D(:,3);
PO2_60SZA = Data_fig5_1D(:,4);
PO2_80SZA = Data_fig5_1D(:,5);
Data_fig5_LBL = readmatrix('Figure5_LBL.csv');
ZO = Data_fig5_LBL(1:30,2);
Press_LBL = interp1(alt,Press_1PAL,ZO);
PO2_0SZA_LBL = Data_fig5_LBL(1:30,3);
Z1 = Data_fig5_LBL(1:30,4);
Press1_LBL = interp1(alt,Press_1PAL,Z1);
PO2_60SZA_LBL = Data_fig5_LBL(1:30,5);
Z2 = Data_fig5_LBL(:,6);
Press2_LBL = interp1(alt,Press_1PAL,Z2);
PO2_80SZA_LBL = Data_fig5_LBL(:,7);
%% Figure 5. COMPARISON WITH FERNANDEZ ET AL., 2006 LINE BY LINE CALCULATION
h1=zeros(5,1);
figure
h1(1)=plot(PO2_0SZA,Press_1D,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
hold on
h1(2)=plot(PO2_60SZA,Press_1D,'Color',[0 0.4470 0.7410],'LineWidth',2);
h1(3)=plot(PO2_80SZA,Press_1D,'Color',[0.8500 0.3250 0.0980],'LineWidth',2);
h1(4)=plot(PO2_0SZA,Press_1D,'Color',[0.4660 0.6740 0.1880],'LineWidth',2);
h1(5)=plot(PO2_0SZA_LBL,Press_LBL,'-.','Color',[0.4660 0.6740 0.1880],'LineWidth',2);

plot(PO2_0SZA,Press_1D,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(PO2_60SZA,Press_1D,'Color',[0 0.4470 0.7410],'LineWidth',2)
plot(PO2_80SZA,Press_1D,'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(PO2_0SZA_LBL,Press_LBL,'-.','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(PO2_60SZA_LBL,Press1_LBL,'-.','Color',[0 0.4470 0.7410],'LineWidth',2)
plot(PO2_80SZA_LBL,Press2_LBL,'-.','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
xlabel('Photolysis rates of O_2 (s^{-1})','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
legend(h1,'SZA = 0','SZA = 60^{\circ}','SZA = 80^{\circ}','Our model','Fernandez et al., 2006','fontsize',14);
set(gca,'XScale','log');set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
set(gca,'Ydir','reverse');
xlim([1e-13,1e-6]);ylim([1E-3 10]);
hold off
%% Figure 6. JO2 AND JH2O comparisons between the 1-D and WACCM3D models
hg_c=zeros(6,1);
figure
subplot(121)
hg_c(1) = plot(PO2_1PAL,Press_1PAL,'Color',[0.47,0.67,0.19],'LineWidth',2);
hold on
hg_c(2) = plot(PO2_01PAL,Press_01PAL,'Color',[0.85,0.33,0.10],'LineWidth',2);
hg_c(3) = plot(PO2_001PAL,Press_001PAL,'Color',[0.93,0.69,0.13],'LineWidth',2);
hg_c(4) = plot(PO2_0001PAL,Press_0001PAL,'Color',[0.00,0.45,0.74],'LineWidth',2);
hg_c(5) = plot(PO2_1PAL,Press_1PAL,'Color',[0.47,0.67,0.19],'LineWidth',2);
hg_c(6) = plot(PO2_1PAL_WACCM3D,TotalP_1PAL,'-.','Color',[0.47,0.67,0.19],'LineWidth',2);

plot(PO2_1PAL,Press_1PAL,'Color',[0.47,0.67,0.19],'LineWidth',2)
hold on
plot(PO2_01PAL,Press_01PAL,'Color',[0.85,0.33,0.10],'LineWidth',2)
plot(PO2_001PAL,Press_001PAL,'Color',[0.93,0.69,0.13],'LineWidth',2)
plot(PO2_0001PAL,Press_0001PAL,'Color',[0.00,0.45,0.74],'LineWidth',2)

plot(PO2_1PAL_WACCM3D,TotalP_1PAL,'-.','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(PO2_01PAL_WACCM3D,TotalP_01PAL,'-.','Color',[0.85,0.33,0.10],'LineWidth',2)
plot(PO2_001PAL_WACCM3D,TotalP_001PAL,'-.','Color',[0.93,0.69,0.13],'LineWidth',2)
plot(PO2_0001PAL_WACCM3D,TotalP_0001PAL,'-.','Color',[0.00,0.45,0.74],'LineWidth',2)
xlim([1e-20,1e-5]);
ylim([1e-2 1e3]);

xlabel('Photolysis requencies of O_2 (s^{-1})','fontweight','bold');
ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'Ydir','reverse');
set(gca,'LineWidth',2,'FontSize',20);
hold off

subplot(122)
plot(PH2O_1PAL,Press_1PAL,'Color',[0.47,0.67,0.19],'LineWidth',2)
hold on
plot(PH2O_01PAL,Press_01PAL,'Color',[0.85,0.33,0.10],'LineWidth',2)
plot(PH2O_001PAL,Press_001PAL,'Color',[0.93,0.69,0.13],'LineWidth',2)
plot(PH2O_0001PAL,Press_0001PAL,'Color',[0.00,0.45,0.74],'LineWidth',2)

plot(PH2O_1PAL_WACCM3D,TotalP_1PAL,'-.','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(PH2O_01PAL_WACCM3D,TotalP_01PAL,'-.','Color',[0.85,0.33,0.10],'LineWidth',2)
plot(PH2O_001PAL_WACCM3D,TotalP_001PAL,'-.','Color',[0.93,0.69,0.13],'LineWidth',2)
plot(PH2O_0001PAL_WACCM3D,TotalP_0001PAL,'-.','Color',[0.00,0.45,0.74],'LineWidth',2)
xlim([1e-20,1e-5]);
ylim([1e-2 1e3]);

xlabel('Photolysis frequencies of H_2O (s^{-1})','fontweight','bold');
ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'Ydir','reverse');
set(gca,'LineWidth',2,'FontSize',20);
hold off
legend(hg_c,'1PAL','0.1PAL','0.01PAL','0.001PAL','1D Gauss','3D from Cooke et al.','fontsize',14)

%% Figure 7: Compare Ox production (cm^-3 s^-1) in the 1D and WACCM3D models
figure
subplot(221)
plot(2.*ProOx_1PAL_1D./1e5,Press_1PAL,'Color',[0.47,0.67,0.19],'LineWidth',2)
hold on
plot(2.*ProOx_NOS_1PAL./1e5,Press_1PAL,'-.','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(2.*ProOx_1PAL_WACCM3D./1e5,TotalP_1PAL,'Color',[0.00,0.45,0.74],'LineWidth',2)
plot(ProOx_1PAL_Kinetics./1e5,Press_Kinetics,'Color',[0.85,0.33,0.10],'LineWidth',2)

ylim([1 1e3]);
set(gca, 'YDir','reverse');
xlabel('Production rate of Ox (10^5 cm^{-3} s^{-1}) at 1 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'YScale','log');
legend('1D','1D only without scattering','WACCM6 3D','Kinetics 1D model');
set(gca,'LineWidth',2,'FontSize',20);
hold off

subplot(222)
plot(2.*ProOx_01PAL_1D./1e5,Press_01PAL,'Color',[0.47,0.67,0.19],'LineWidth',2)
hold on
plot(2.*ProOx_NOS_01PAL./1e5,Press_01PAL,'-.','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(2.*ProOx_01PAL_WACCM3D./1e5,TotalP_01PAL,'Color',[0.00,0.45,0.74],'LineWidth',2)
plot(ProOx_01PAL_Kinetics./1e5,Press_Kinetics,'Color',[0.85,0.33,0.10],'LineWidth',2)

ylim([1 1e3]);
set(gca, 'YDir','reverse');
xlabel('Production rate of Ox (10^5 cm^{-3} s^{-1}) at 0.1 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off

subplot(223)
plot(2.*ProOx_001PAL_1D./1e5,Press_001PAL,'Color',[0.47,0.67,0.19],'LineWidth',2)
hold on
plot(2.*ProOx_NOS_001PAL./1e5,Press_001PAL,'-.','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(2.*ProOx_001PAL_WACCM3D./1e5,TotalP_001PAL,'Color',[0.00,0.45,0.74],'LineWidth',2);
plot(ProOx_001PAL_Kinetics./1e5,Press_Kinetics,'Color',[0.85,0.33,0.10],'LineWidth',2)

ylim([1 1e3]);
set(gca, 'YDir','reverse');
xlabel('Production rate of Ox (10^5 cm^{-3} s^{-1}) at 0.01 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off

subplot(224)
plot(2.*ProOx_0001PAL_1D./1e5,Press_0001PAL,'Color',[0.47,0.67,0.19],'LineWidth',2);
hold on
plot(2.*ProOx_NOS_0001PAL./1e5,Press_0001PAL,'-.','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(2.*ProOx_0001PAL_WACCM3D./1e5,TotalP_0001PAL,'Color',[0.00,0.45,0.74],'LineWidth',2);
plot(ProOx_0001PAL_Kinetics./1e5,Press_Kinetics,'Color',[0.85,0.33,0.10],'LineWidth',2)

ylim([1 1e3]);
set(gca, 'YDir','reverse');
xlabel('Production rate of Ox (10^5 cm^{-3} s^{-1}) at 0.001 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off

%% 1.6 Compare JO2 & JH2O with/without scattering in the 1D
% *************** READ DATA FROM 'O3columndepth' **************
Data_fig8_1D = readmatrix('Figure8_1D');
% Pressure (hPa)
Press_fig8_1D = Data_fig8_1D(:,2); 
% photolysis rates with scattering
PO2_1PAL_S = Data_fig8_1D(:,3); 
PO2_001PAL_S = Data_fig8_1D(:,4); 
PO2_00001PAL_S = Data_fig8_1D(:,5); 
PH2O_1PAL_S = Data_fig8_1D(:,6); 
PH2O_001PAL_S = Data_fig8_1D(:,7); 
PH2O_00001PAL_S = Data_fig8_1D(:,8); 
%% READ PHOTOLYSIS RATES FOR O2 AND H2O ONLY ONE STEP WITHOUT TWO STREAM
PO2_1PAL_NO = Data_fig8_1D(:,9); 
PO2_001PAL_NO = Data_fig8_1D(:,10); 
PO2_00001PAL_NO = Data_fig8_1D(:,11); 
PH2O_1PAL_NO = Data_fig8_1D(:,12); 
PH2O_001PAL_NO = Data_fig8_1D(:,13); 
PH2O_00001PAL_NO = Data_fig8_1D(:,14); 
% WITH BAND MODEL
PO2_1PAL_BAND = Data_fig8_1D(:,15); 
PO2_001PAL_BAND = Data_fig8_1D(:,16); 
PO2_00001PAL_BAND = Data_fig8_1D(:,17); 
PH2O_1PAL_BAND = Data_fig8_1D(:,18); 
PH2O_001PAL_BAND = Data_fig8_1D(:,19); 
PH2O_00001PAL_BAND = Data_fig8_1D(:,20); 

%% Figure 8 PLOT PO2 AT 1PAL WITH/WITHOUT SCATTERING/BAND in pressure
figure
subplot(121)
plot(PO2_1PAL,Press_1PAL,'Color',[0.4940 0.1840 0.5560],'LineWidth',2)
hold on
plot(PO2_1PAL_NO,Press_1D,'-.','Color',[0.4940 0.1840 0.5560],'LineWidth',2)
plot(PO2_1PAL_BAND,Press_1D,':','Color',[0.4940 0.1840 0.5560],'LineWidth',2)
plot(PO2_001PAL_S,Press_1D,'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(PO2_001PAL_NO,Press_1D,'-.','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(PO2_00001PAL_S,Press_1D,'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
plot(PO2_00001PAL_NO,Press_1D,'-.','Color',[0.9290 0.6940 0.1250],'LineWidth',2)
xlim([1e-15 1e-7]);ylim([1e-3 1e3]);
xlabel('Photolysis rates of O_2 (s^{-1})','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
legend('1 PAL with scattering','1 PAL without scattering','1 PAL with Band model','10^{-2} PAL with scattering','10^{-2} PAL without scattering','10^{-4} PAL with scattering','10^{-4} PAL without scattering','FontSize',14);
set(gca,'XScale','log');set(gca,'YScale','log');
set(gca, 'YDir','reverse');
set(gca,'LineWidth',2,'FontSize',20);
hold off

subplot(122)
plot(PH2O_1PAL_S,Press_1D,'Color',[0.4940 0.1840 0.5560],'LineWidth',2)
hold on
plot(PH2O_1PAL_NO,Press_1D,'-.','Color',[0.4940 0.1840 0.5560],'LineWidth',2)
plot(PH2O_1PAL_BAND,Press_1D,':','Color',[0.4940 0.1840 0.5560],'LineWidth',2)
plot(PH2O_001PAL_S,Press_1D,'Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(PH2O_001PAL_NO,Press_1D,'-.','Color',[0.8500 0.3250 0.0980],'LineWidth',2)
plot(PH2O_00001PAL_S,Press_1D,'Color',[0.9290 0.6940 0.1250],'LineWidth',2)
plot(PH2O_00001PAL_NO,Press_1D,'-.','Color',[0.9290 0.6940 0.1250],'LineWidth',2)
xlim([1e-15 1e-5]);ylim([1e-3 1e3]);
xlabel('Photolysis rates of H_2O (s^{-1})','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
legend('1 PAL with scattering','1 PAL without scattering','1 PAL with band model','10^{-2} PAL with scattering','10^{-2} PAL without scattering','10^{-4} PAL with scattering','10^{-4} PAL without scattering','FontSize',14);
set(gca,'XScale','log');set(gca,'YScale','log');
set(gca, 'YDir','reverse');
set(gca,'LineWidth',2,'FontSize',20);
hold off

%% Figure 9: compare NumOH at 1 PAL, 0.1 PAL, 0.01 PAL and 0.001 PAL
figure
subplot(221)
plot(NumOH_1PAL,Press_1PAL,'Color',[0.47,0.67,0.19],'LineWidth',2)
hold on
plot(NOH_1PAL_WACCM3D,NP_1PAL_WACCM3D,'Color',[0.00,0.45,0.74],'LineWidth',2);
plot(NumOH_ROCK3D(:,1),Press_ROCK3D,'Color',[0.93,0.69,0.13],'LineWidth',2)
xlim([6e5 2e6]);
ylim([200 950]);
set(gca, 'YDir','reverse');
xlabel('OH number density (cm^{-3}) at 1 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'XScale','log');
set(gca,'LineWidth',2,'FontSize',20);
legend('Standard 1D','WACCM6 3D','ROCKE 3D','fontsize',14);
hold off

subplot(222)
plot(NumOH_01PAL,Press_01PAL,'Color',[0.47,0.67,0.19],'LineWidth',2)
hold on
plot(NOH_01PAL_WACCM3D,NP_01PAL_WACCM3D,'Color',[0.00,0.45,0.74],'LineWidth',2);
plot(NumOH_ROCK3D(:,2),Press_ROCK3D,'Color',[0.93,0.69,0.13],'LineWidth',2)
xlim([2e5 1e7]);
ylim([200 950]);
set(gca, 'YDir','reverse');
xlabel('OH number density (cm^{-3}) at 0.1 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'XScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off

subplot(223)
plot(NumOH_001PAL,Press_001PAL,'Color',[0.47,0.67,0.19],'LineWidth',2)
hold on
plot(NOH_001PAL_WACCM3D,NP_1PAL_WACCM3D,'Color',[0.00,0.45,0.74],'LineWidth',2);
plot(NumOH_ROCK3D(:,3),Press_ROCK3D,'Color',[0.93,0.69,0.13],'LineWidth',2)
xlim([2e5 3e8]);
ylim([200 950]);
set(gca, 'YDir','reverse');
xlabel('OH number density (cm^{-3}) at 0.01 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'XScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off

subplot(224)
plot(NumOH_0001PAL,Press_0001PAL,'Color',[0.47,0.67,0.19],'LineWidth',2)
hold on
plot(NOH_0001PAL_WACCM3D,NP_0001PAL_WACCM3D,'Color',[0.00,0.45,0.74],'LineWidth',2);
plot(NumOH_ROCK3D(:,4),Press_ROCK3D,'Color',[0.93,0.69,0.13],'LineWidth',2)
xlim([2e5 3e8]);
ylim([200 950]);
set(gca, 'YDir','reverse');
xlabel('OH number density (cm^{-3}) at 0.001 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'XScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off

%% Figure 10 PLOT REACTION RATES OF OH
figure
subplot(221)
plot(R1OH_1PAL,Press_1PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(R25OH_1PAL,Press_1PAL,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(ProOH_1PAL_3D,TotalP_1PAL_3D,'Color',[0.00,0.45,0.74],'LineWidth',2)
plot(ROH1PAL_WACCM3D,TotalP_1PAL,'--','Color',[0.00,0.45,0.74],'LineWidth',2)
xlabel('Reaction rate at 1 PAL (cm^{-3} s^{-1})','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
legend('O(^1D) + H_2O --> 2OH in 1D','H_2O + hv --> H + OH in 1D','O(^1D) + H_2O --> 2OH in WACCM 3D','H_2O + hv --> H + OH in WACCM 3D','Fontsize',10);
xlim([1 1e9]);
ylim([1e-2 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

subplot(222)
plot(R1OH_01PAL,Press_01PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(R25OH_01PAL,Press_01PAL,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(ProOH_01PAL_3D,TotalP_01PAL_3D,'Color',[0.00,0.45,0.74],'LineWidth',2)
plot(ROH01PAL_WACCM3D,TotalP_01PAL,'--','Color',[0.00,0.45,0.74],'LineWidth',2)
xlabel('Reaction rate at 0.1 PAL (cm^{-3} s^{-1})','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
xlim([1 1e9]);
ylim([1e-2 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

subplot(223)
plot(R1OH_001PAL,Press_001PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(R25OH_001PAL,Press_001PAL,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(ProOH_001PAL_3D,TotalP_001PAL_3D,'Color',[0.00,0.45,0.74],'LineWidth',2)
plot(ROH001PAL_WACCM3D,TotalP_001PAL,'--','Color',[0.00,0.45,0.74],'LineWidth',2)
xlabel('Reaction rate at 0.01 PAL (cm^{-3} s^{-1})','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
xlim([1 1e9]);
ylim([1e-2 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

subplot(224)
plot(R1OH_0001PAL,Press_0001PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(R25OH_0001PAL,Press_0001PAL,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(ProOH_0001PAL_3D,TotalP_0001PAL_3D,'Color',[0.00,0.45,0.74],'LineWidth',2)
plot(ROH0001PAL_WACCM3D,TotalP_0001PAL,'--','Color',[0.00,0.45,0.74],'LineWidth',2)

xlabel('Reaction rate at 0.001 PAL (cm^{-3} s^{-1})','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
xlim([1 1e9]);
ylim([1e-2 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

%% Figure 11. PLOT REACTION RATES OF OH from H2O photolysis in 1D and WACCM 3D
figure
subplot(221)
plot(R25OH_1PAL,Press_1PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(R25OH_NO_1PAL,Press_1PAL,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(ROH1PAL_WACCM3D,TotalP_1PAL,'--','Color',[0.00,0.45,0.74],'LineWidth',2)

xlabel('Reaction rate at 1 PAL (cm^{-3} s^{-1})','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
legend('H_2O + hv --> H + OH in 1D','H_2O + hv --> H + OH in 1D no S & A','H_2O + hv --> H + OH in WACCM 3D','Fontsize',14);
xlim([1 1e9]);
ylim([1e-2 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

subplot(222)
plot(R25OH_01PAL,Press_01PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(R25OH_NO_01PAL,Press_01PAL,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(ROH01PAL_WACCM3D,TotalP_01PAL,'--','Color',[0.00,0.45,0.74],'LineWidth',2)
plot(RateOH_01PAL_WACCM1206,Press_01PAL_WACCM1206,'Color',[0.00,0.45,0.74],'LineWidth',2)

xlabel('Reaction rate at 0.1 PAL (cm^{-3} s^{-1})','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
xlim([1 1e9]);
ylim([1e-2 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

subplot(223)
plot(R25OH_001PAL,Press_001PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(R25OH_NO_001PAL,Press_001PAL,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(ROH001PAL_WACCM3D,TotalP_001PAL,'--','Color',[0.00,0.45,0.74],'LineWidth',2)
plot(RateOH_001PAL_WACCM1206,Press_001PAL_WACCM1206,'Color',[0.00,0.45,0.74],'LineWidth',2)

xlabel('Reaction rate at 0.01 PAL (cm^{-3} s^{-1})','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
xlim([1 1e9]);
ylim([1e-2 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

subplot(224)
plot(R25OH_0001PAL,Press_0001PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(R25OH_NO_0001PAL,Press_0001PAL,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(ROH0001PAL_WACCM3D,TotalP_0001PAL,'--','Color',[0.00,0.45,0.74],'LineWidth',2)
plot(RateOH_0001PAL_WACCM1206,Press_0001PAL_WACCM1206,'Color',[0.00,0.45,0.74],'LineWidth',2)

xlabel('Reaction rate at 0.001 PAL (cm^{-3} s^{-1})','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
xlim([1 1e9]);
ylim([1e-2 1e3]);
legend('H_2O + hv --> H + OH in 1D','H_2O + hv --> H + OH in 1D no S & A','H_2O + hv --> H + OH in WACCM 3D','Fontsize',14);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

%% Figure 12 OH number density with/without SR in the 1D and WACCM6
figure
subplot(221)
plot(NumOH_1PAL,Press_1PAL,'Color',[0.47,0.67,0.19],'LineWidth',2)
hold on
plot(NumOH_NOS_1PAL,Press_1PAL,'-.','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(NumOH_NO_1PAL,Press_1PAL,'--','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(NOH_1PAL_WACCM3D,NP_1PAL_WACCM3D,'--','Color',[0.00,0.45,0.74],'LineWidth',2);
set(gca, 'YDir','reverse');
%xlim([1e5 4.5e8]);
ylim([200 950]);
xlabel('OH number density (cm^{-3}) at 1 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
legend('Standard 1D','No scattering 1D','No scattering and absorbers 1D','WACCM6','fontsize',14);
set(gca,'XScale','log');%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off
subplot(222)
plot(NumOH_01PAL,Press_01PAL,'Color',[0.47,0.67,0.19],'LineWidth',2)
hold on
plot(NumOH_NOS_01PAL,Press_01PAL,'-.','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(NumOH_NO_01PAL,Press_01PAL,'--','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(NOH_01PAL_WACCM3D,NP_01PAL_WACCM3D,'--','Color',[0.00,0.45,0.74],'LineWidth',2);
plot(NumOH_01PAL_WACCM1206,Press_01PAL_WACCM1206,'Color',[0.00,0.45,0.74],'LineWidth',2);
set(gca, 'YDir','reverse');
xlim([2e5 1e7]);
ylim([200 950]);
xlabel('OH number density (cm^{-3}) at 0.1 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'XScale','log');%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off

subplot(223)
plot(NumOH_001PAL,Press_001PAL,'Color',[0.47,0.67,0.19],'LineWidth',2)
hold on
plot(NumOH_NOS_001PAL,Press_001PAL,'-.','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(NumOH_NO_001PAL,Press_001PAL,'--','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(NOH_001PAL_WACCM3D,NP_001PAL_WACCM3D,'--','Color',[0.00,0.45,0.74],'LineWidth',2);
plot(NumOH_001PAL_WACCM1206,Press_001PAL_WACCM1206,'Color',[0.00,0.45,0.74],'LineWidth',2);

set(gca, 'YDir','reverse');
xlim([2e5 3e8]);
ylim([200 950]);
xlabel('OH number density (cm^{-3}) at 0.01 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'XScale','log');%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off

subplot(224)
plot(NumOH_0001PAL,Press_0001PAL,'Color',[0.47,0.67,0.19],'LineWidth',2)
hold on
plot(NumOH_NOS_0001PAL,Press_0001PAL,'-.','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(NumOH_NO_0001PAL,Press_0001PAL,'--','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(NOH_0001PAL_WACCM3D,NP_0001PAL_WACCM3D,'--','Color',[0.00,0.45,0.74],'LineWidth',2);
plot(NumOH_0001PAL_WACCM1206,Press_0001PAL_WACCM1206,'Color',[0.00,0.45,0.74],'LineWidth',2);

set(gca, 'YDir','reverse');
xlim([2e5 4.5e8]);
ylim([200 950]);
xlabel('OH number density (cm^{-3}) at 0.001 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'XScale','log');%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off
%% Data for plotting surface lifetime vs. column integrated lifetime for each CH4 flux

% Present Methane FLUX
lifetime_surface = [4.6 9.14 6.16 2.10 0.745];
%lifetime = [7.28, 12.84 7.05 1.46 0.54];
lifetime = [8.1 14.3 7.69 1.95 0.68];
lifetime_cooke = [9 5 0.6 0.1]; 
surf_CH4_1D = [1.6 2.68 1.3 0.28 0.1]*1e-6;
pO2_WACCM = [0.1 0.01];
surf_CH4_WACCM = [1.2 7.6E-2]*1e-6;

%% Figure 13 plot surface lifetime vs. column integrated lifetime for each CH4 flux
figure
subplot(121)
plot(pO2,surf_CH4_1D,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(pO2(1:4),[0.8 0.8 0.8 0.8]*1e-6,'Color',[0.93,0.69,0.13],'LineWidth',2)
plot(pO2_WACCM,surf_CH4_WACCM,'--','Color',[0.93,0.69,0.13],'LineWidth',2)
xlabel('pO_2 (PAL)','fontweight','bold');
ylabel('Surface CH_4 mixing ratio','fontweight','bold');
ylim([5e-8 3e-6]);
legend('1-D','Fixed CH_4 mixing ratio in WACCM6','Fixed CH_4 flux in WACCM6','Fontsize',10);
set(gca,'XScale','log');set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);

subplot(122)
plot(pO2,lifetime,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(pO2(1:4),lifetime_cooke,'Color',[0.93,0.69,0.13],'LineWidth',2)
xlabel('pO_2 (PAL)','fontweight','bold');
ylabel('CH_4 lifetime (yr)','fontweight','bold');
set(gca,'XScale','log');set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
%% SUPP. INFO.
% ****** plot: O3 column depth vs. pO2 for different calculations in the 1D******
figure
plot(pO2,Gau_LOWH2O_O3,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(pO2,Gau_O3,'--*','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(pO2,new_O3,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(pO2,old_O3,'-.','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(pO2,O3_nolit,'--o','Color',[0.4660 0.6740 0.1880],'LineStyle',':','LineWidth',2)
plot(pO2(1:4),Even_O3,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(pO2,ones(5)*292,':','LineWidth',2)
xlabel('pO_2 (PAL)','fontweight','bold');ylabel('O_3 column depth (DU)','fontweight','bold');
text(1e-3,340,'WOUDC Average = 292 DU','fontsize',14);
legend('8 Gauss points','8 Gauss points with higher H_2O','ZY = 48.2^{\circ}','ZY = 48.2^{\circ}, Old H_2O cross-sections','ZY = 48.2^{\circ}, No lightning','Even distributed','fontsize',14,'Location','best');
set(gca,'XScale','log');%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
ylim([0 350]);
%% DENSITY PROFILE
% ***** plot figure S2. ***************
str = zeros(length(TotalP_1PAL),1);
figure
plot((DEN_1PAL_WACCM3D1-DEN_1PAL_WACCM3D)./DEN_1PAL_WACCM3D*100,TotalP_1PAL,'Color',[0.47,0.67,0.19],'LineWidth',2)
hold on
%plot(DEN_1PAL_WACCM3D1,TotalP_1PAL,'--','Color',[0.00,0.45,0.74],'LineWidth',2)
plot((DEN_01PAL_WACCM3D1-DEN_01PAL_WACCM3D)./DEN_01PAL_WACCM3D*100,TotalP_01PAL,'Color',[0.85,0.33,0.10],'LineWidth',2)
%plot(DEN_01PAL_WACCM3D1,TotalP_01PAL,'--','Color',[0.93,0.69,0.13],'LineWidth',2)
plot((DEN_001PAL_WACCM3D1-DEN_001PAL_WACCM3D)./DEN_001PAL_WACCM3D*100,TotalP_001PAL,'Color',[0.93,0.69,0.13],'LineWidth',2)
%plot(DEN_001PAL_WACCM3D1,TotalP_001PAL,'--','Color',[0.47,0.67,0.19],'LineWidth',2)
plot((DEN_0001PAL_WACCM3D1-DEN_0001PAL_WACCM3D)./DEN_0001PAL_WACCM3D*100,TotalP_0001PAL,'Color',[0.00,0.45,0.74],'LineWidth',2)
%plot(DEN_0001PAL_WACCM3D1,TotalP_0001PAL,'--','Color',[0.85,0.33,0.10],'LineWidth',2)
plot(str,TotalP_1PAL,'--','Color',[0,0,0],'LineWidth',2)

%xlabel('Total number density (cm^{-3})');
xlabel('Difference to the density in WACCM 3D (%)');
ylabel('Pressure (hPa)');
ylim([1e-2 1e3]);
legend('1 PAL','0.1 PAL','0.01 PAL','0.001 PAL');
set(gca, 'YDir','reverse');
set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
%%
PO2_1PAL = xlsread('O3columnDepth','48_2Degree','B2:B35');
PO2_01PAL = xlsread('O3columnDepth','48_2Degree','C2:C35');
PO2_001PAL = xlsread('O3columnDepth','48_2Degree','D2:D35');
PO2_0001PAL = xlsread('O3columnDepth','48_2Degree','E2:E35');
PO2_00001PAL = xlsread('O3columnDepth','48_2Degree','F2:F35');
PH2O_1PAL = xlsread('O3columnDepth','48_2Degree','G2:G35');
PH2O_01PAL = xlsread('O3columnDepth','48_2Degree','H2:H35');
PH2O_001PAL = xlsread('O3columnDepth','48_2Degree','I2:I35');
PH2O_0001PAL = xlsread('O3columnDepth','48_2Degree','J2:J35');
PH2O_00001PAL = xlsread('O3columnDepth','48_2Degree','K2:K35');
% with 60 degree solar zenith angle
SPO2_1PAL = xlsread('O3columnDepth','60_Degree','B2:B35');
SPO2_01PAL = xlsread('O3columnDepth','60_Degree','C2:C35');
SPO2_001PAL = xlsread('O3columnDepth','60_Degree','D2:D35');
SPO2_0001PAL = xlsread('O3columnDepth','60_Degree','E2:E35');
SPO2_00001PAL = xlsread('O3columnDepth','60_Degree','F2:F35');
SPH2O_1PAL = xlsread('O3columnDepth','60_Degree','G2:G35');
SPH2O_01PAL = xlsread('O3columnDepth','60_Degree','H2:H35');
SPH2O_001PAL = xlsread('O3columnDepth','60_Degree','I2:I35');
SPH2O_0001PAL = xlsread('O3columnDepth','60_Degree','J2:J35');
SPH2O_00001PAL = xlsread('O3columnDepth','60_Degree','K2:K35');
% with 8 gauss points
PZ = xlsread('O3columnDepth','GaussPoints','A2:A35');
GPO2_1PAL = xlsread('O3columnDepth','GaussPoints','B2:B35');
GPO2_01PAL = xlsread('O3columnDepth','GaussPoints','C2:C35');
GPO2_001PAL = xlsread('O3columnDepth','GaussPoints','D2:D35');
GPO2_0001PAL = xlsread('O3columnDepth','GaussPoints','E2:E35');
GPO2_00001PAL = xlsread('O3columnDepth','GaussPoints','F2:F35');
GPH2O_1PAL = xlsread('O3columnDepth','GaussPoints','G2:G35');
GPH2O_01PAL = xlsread('O3columnDepth','GaussPoints','H2:H35');
GPH2O_001PAL = xlsread('O3columnDepth','GaussPoints','I2:I35');
GPH2O_0001PAL = xlsread('O3columnDepth','GaussPoints','J2:J35');
GPH2O_00001PAL = xlsread('O3columnDepth','GaussPoints','K2:K35');
% **************** plot figure S3 *********************
% compare 48.2, 60 degree and 8 gauss poionts at 1 0.01 0.0001PAL
h_gauss=zeros(6,1);
figure
subplot(211)
h_gauss(1)=plot(GPO2_1PAL,PZ./1e5,'Color',[0.00,0.45,0.74],'LineWidth',2);
hold on
h_gauss(2)=plot(PO2_1PAL,PZ./1e5,'--','Color',[0.00,0.45,0.74],'LineWidth',2);
h_gauss(3)=plot(SPO2_1PAL,PZ./1e5,'-.','Color',[0.00,0.45,0.74],'LineWidth',2);
h_gauss(4)=plot(SPO2_1PAL,PZ./1e5,'-.','Color',[0.00,0.45,0.74],'LineWidth',2);
h_gauss(5)=plot(SPO2_001PAL,PZ./1e5,'-.','Color',[0.93,0.69,0.13],'LineWidth',2);
h_gauss(6)=plot(SPO2_00001PAL,PZ./1e5,'-.','Color',[0.47,0.67,0.19],'LineWidth',2);

plot(GPO2_1PAL,PZ./1e5,'Color',[0.00,0.45,0.74],'LineWidth',2)
hold on
plot(GPO2_001PAL,PZ./1e5,'Color',[0.93,0.69,0.13],'LineWidth',2)
plot(GPO2_00001PAL,PZ./1e5,'Color',[0.47,0.67,0.19],'LineWidth',2)
plot(PO2_1PAL,PZ./1e5,'--','Color',[0.00,0.45,0.74],'LineWidth',2)
plot(PO2_001PAL,PZ./1e5,'--','Color',[0.93,0.69,0.13],'LineWidth',2)
plot(PO2_00001PAL,PZ./1e5,'--','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(SPO2_1PAL,PZ./1e5,'-.','Color',[0.00,0.45,0.74],'LineWidth',2)
plot(SPO2_001PAL,PZ./1e5,'-.','Color',[0.93,0.69,0.13],'LineWidth',2)
plot(SPO2_00001PAL,PZ./1e5,'-.','Color',[0.47,0.67,0.19],'LineWidth',2)
xlim([1e-20,1e-7]);ylim([0 100]);
xlabel('Photolysis rates of O_2 (s^{-1})','fontweight','bold');ylabel('Altitude (km)','fontweight','bold');
set(gca,'XScale','log');%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off

subplot(212)
plot(GPH2O_1PAL,PZ./1e5,'Color',[0.00,0.45,0.74],'LineWidth',2)
hold on
plot(GPH2O_001PAL,PZ./1e5,'Color',[0.93,0.69,0.13],'LineWidth',2)
plot(GPH2O_00001PAL,PZ./1e5,'Color',[0.47,0.67,0.19],'LineWidth',2)
plot(PH2O_1PAL,PZ./1e5,'--','Color',[0.00,0.45,0.74],'LineWidth',2)
plot(PH2O_001PAL,PZ./1e5,'--','Color',[0.93,0.69,0.13],'LineWidth',2)
plot(PH2O_00001PAL,PZ./1e5,'--','Color',[0.47,0.67,0.19],'LineWidth',2)
plot(SPH2O_1PAL,PZ./1e5,'-.','Color',[0.00,0.45,0.74],'LineWidth',2)
plot(SPH2O_001PAL,PZ./1e5,'-.','Color',[0.93,0.69,0.13],'LineWidth',2)
plot(SPH2O_00001PAL,PZ./1e5,'-.','Color',[0.47,0.67,0.19],'LineWidth',2)
xlim([1e-20,1e-5]);
xlabel('Photolysis rates of H_2O (s^{-1})','fontweight','bold');ylabel('Altitude (km)','fontweight','bold');
set(gca,'XScale','log');%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off
legend(h_gauss,'8 Gauss points','ZY=48.2^{\circ}','ZY=60^{\circ}','1PAL','0.01PAL','0.0001PAL','fontsize',14)

%% global integration with Gaussian Quadrature
% ************ plot figure S5.
Ngau=[2,3,4,5,6,8,10,12];
GPO2_M=[1.50E-9,1.55E-9,1.56E-9,1.57E-9,1.58E-9,1.58E-9,1.59E-9,1.59E-9];
GPO2_T=[5.75E-8,5.94E-8,6.01E-8,6.04E-8,6.06E-8,6.08E-8,6.09E-8,6.09E-8];
%ZY=[48.2,60];AGL=[0.375,0.5];
PO2_M=[5.11E-10,6.17E-10];
PO2_T=[1.92E-8,2.21E-8];
% plot
figure
subplot(211)
plot(Ngau,GPO2_M,'LineWidth',2)
hold on
plot(Ngau,ones(8)*1.59E-9,'-.','LineWidth',2)
text(6,1.55e-9,'Altitude = 60 km','FontSize',20);
%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off

subplot(212)
plot(Ngau,GPO2_T,'LineWidth',2)
hold on
plot(Ngau,ones(8)*6.09E-8,'-.','LineWidth',2)
%ylim([1e-9,6.5e-8]);
text(6,6.e-8,'Altitude = 99.5 km','FontSize',20);
xlabel('Number of Gauss points','fontweight','bold');
ylabel('Photolysis rates of O_2 (s^{-1})','fontweight','bold');
%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off
%% PRODUCTION OF NO FROM LIGHTNING
% ************ plot Figure S6.
PNONOW = 3.574E-2;
PNO = [3.565E-2 9.853E-3 3.67E-3 2.886E-3 2.805E-3];
Prod_NO = 3.E9*(PNO./PNONOW);
figure
plot(pO2,Prod_NO,'LineWidth',2)
xlabel('pO_2 (PAL)','Fontweight','bold');ylabel('Production rate of NO from lightning (cm^{-2} s^{-1})','Fontweight','bold');
set(gca,'XScale','log');set(gca,'YScale','log');ylim([1e8 1e10]);
set(gca,'LineWidth',2,'FontSize',20);


%% Compare solar flux
% read solar flux in the 1-D model
EFLUX = readmatrix('solarflux_1d_1PAL.csv');
WL1 = EFLUX(1:10,2); %Wavelength for farUV
WL2 = EFLUX(11:118,2); %Wavelength for short wavelength and long wavelength
EFLUX_WL1 = EFLUX(1:10,3);
EFLUX_WL2 = EFLUX(11:118,3);
% read solar flux in the 1-D model with farUV from WACCM
EFLUX0726 = readmatrix('SFX_1d_far3d_1PAL0726.csv');
WL1_0726 = EFLUX0726(1:10,2); %Wavelength for farUV
WL2_0726 = EFLUX0726(11:118,2); %Wavelength for short wavelength and long wavelength
EFLUX_WL1_0726 = EFLUX0726(1:10,3);
EFLUX_WL2_0726 = EFLUX0726(11:118,3);
%
EFLUX0726NEW = readmatrix('SFX_1d_far3d_1PAL0726_NEW.csv');

%% get the plots
figure
subplot(2,2,1)
plot(WL1./10,EFLUX_WL1,'LineWidth',2)
xlabel('Wavelength (nm)','Fontweight','bold');
ylabel('Flux at the top (W/m^2/nm)','Fontweight','bold');
xlim([100 175]);ylim([0 2.e-3]);
%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);

subplot(2,2,2)
plot(WL2./10,EFLUX_WL2,'LineWidth',2)
xlabel('Wavelength (nm)','Fontweight','bold');
ylabel('Flux at the top (W/m^2/nm)','Fontweight','bold');
%set(gca,'YScale','log');
xlim([175 700]);ylim([-0.05 2.5]);
set(gca,'LineWidth',2,'FontSize',20);

subplot(2,2,3:4)
plot(EFLUX(:,2)./10,EFLUX(:,3),'LineWidth',2)
xlabel('Wavelength (nm)','Fontweight','bold');
ylabel('Flux at the top (W/m^2/nm)','Fontweight','bold');
xlim([100 400]);ylim([-0.05 2.5]);
%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);

%
figure
plot(EFLUX(:,2)./10,EFLUX(:,3)./EFLUX0726(:,3),'LineWidth',2)
xlabel('Wavelength (nm)','Fontweight','bold');
ylabel('Flux ratio (old/new)','Fontweight','bold');
xlim([100 240]);ylim([-0.05 2.5]);
%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);

EFLUX_WACCM = average';
%%
figure
plot(EFLUX(1:10,2)./10,EFLUX0726(1:10,3).*1e3/EFLUX_WACCM,'LineWidth',2)
xlabel('Wavelength (nm)','Fontweight','bold');
ylabel('Flux ratio (old/new)','Fontweight','bold');
xlim([100 240]);ylim([-0.05 2.5]);
%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
%% TRY THE NEW WAY TO COMPARE THE SOLAR FLUX
WAV_WACCM = wavelength_WACCM(121:241);
SFX_WACCM = SFX(121:241);
WAV_1D = EFLUX0726(2:45,2)./10;
SFX_WACCMfor1D = interp1(WAV_WACCM,SFX_WACCM,WAV_1D);
SFXfor1D = EFLUX(2:45,3)*1e3;
ratio_1D_WACCM = SFXfor1D./SFX_WACCMfor1D;

figure
plot(EFLUX(:,2)./10,EFLUX(:,3)./EFLUX0726NEW(:,3),'LineWidth',2)
xlabel('Wavelength (nm)','Fontweight','bold');
ylabel('Flux ratio (old/new)','Fontweight','bold');
xlim([100 240]);
ylim([-0.05 2.5]);
%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);

%% Different Gaussion points VS. 48.2 degree solar zenith angle
Integ_Gauss = readmatrix('GaussianTEST');
figure
plot([8 6 5 4 3 2 1],Integ_Gauss(1:7),'o-','LineWidth',2)
hold on
%plot(1,Integ_Gauss(8),'*', 'MarkerSize', 10)
plot([8 6 5 4 3 2 1],ones(7)*Integ_Gauss(9),'-.','LineWidth',2 )
xlabel('Number of Gaussian Points');
ylabel({'Column-integrated reaction rate of';'O_2 photolysis'});
%text(1.2,Integ_Gauss(8),'\leftarrow 1 Gaussian point with 48.2 degree');
text(4,6.4e12,'Solar zenith angle = 48.2^\circ','FontSize',18);
set(gca,'LineWidth',2,'FontSize',20);
%% Section 2: Discussion about Methane


%% plot for figure S8?
figure
subplot(221)
plot(fCH4_1PAL,press_fH2O_1PAL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
hold on
plot(fCH4_1PAL_WACCM3D,TotalP_1PAL,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(fCH4_ROCK3D(:,1),Press_ROCK3D,'-.','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
xlabel('CH_4 mixing ratio at 1 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
legend('1D','WACCM 3D','ROCKE 3D','Fontsize',14);
xlim([1e-10 1e-5]);
ylim([1e-1 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

subplot(222)
plot(fCH4_01PAL,press_fH2O_01PAL,'Color',[0.85,0.33,0.10],'LineWidth',2)
hold on
plot(fCH4_01PAL_WACCM3D,TotalP_01PAL,'--','Color',[0.85,0.33,0.10],'LineWidth',2)
plot(fCH4_ROCK3D(:,2),Press_ROCK3D,'-.','Color',[0.85,0.33,0.10],'LineWidth',2)
xlabel('CH_4 mixing ratio at 0.1 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
xlim([1e-10 1e-5]);
ylim([1e-1 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

subplot(223)
plot(fCH4_001PAL,press_fH2O_001PAL,'Color',[0.93,0.69,0.13],'LineWidth',2)
hold on
plot(fCH4_001PAL_WACCM3D,TotalP_001PAL,'--','Color',[0.93,0.69,0.13],'LineWidth',2)
plot(fCH4_ROCK3D(:,3),Press_ROCK3D,'-.','Color',[0.93,0.69,0.13],'LineWidth',2)
xlabel('CH_4 mixing ratio at 0.01 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
xlim([1e-10 1e-5]);
ylim([1e-1 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

subplot(224)
plot(fCH4_0001PAL,press_fH2O_0001PAL,'Color',[0.00,0.45,0.74],'LineWidth',2)
hold on
plot(fCH4_0001PAL_WACCM3D,TotalP_0001PAL,'--','Color',[0.00,0.45,0.74],'LineWidth',2)
plot(fCH4_ROCK3D(:,4),Press_ROCK3D,'-.','Color',[0.00,0.45,0.74],'LineWidth',2)
xlabel('CH_4 mixing ratio at 0.0001 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
set(gca,'LineWidth',2,'FontSize',20);
xlim([1e-10 1e-5]);
ylim([1e-1 1e3]);
set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca, 'YDir','reverse');
hold off

%%
figure
plot(PH2O2_01PAL,Press_01PAL,'LineWidth',2)
hold on
plot(PH2O2_NOS_01PAL,Press_01PAL,'-.','LineWidth',2)
plot(PH2O2_NO_01PAL,Press_01PAL,'--','LineWidth',2)
ylim([200 950]);
set(gca, 'YDir','reverse');
xlabel('PH2O2 (cm^{-3} s^{-1}) at 0.1 PAL','fontweight','bold');ylabel('Pressure (hPa)','fontweight','bold');
%set(gca,'XScale','log');
set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
hold off
%% *** plot: O3 column depth vs. pO2 for different calculations in the 1D with CH3CL***
figure
subplot(211)
plot(pO2(1:4),O3_CL,':','Color',[0 0 0],'LineWidth',2)
hold on
plot(pO2(1:4),O3_LOWCL,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
plot(pO2(1:4),O3_CL05ppt,'Color',[0.93,0.69,0.13],'LineWidth',2)
plot(pO2(1:3),O3_CL05pptNOSNOA,'Color',[0.85,0.33,0.10],'LineWidth',2)
plot(pO2,ones(5)*292,':','LineWidth',2)
xlabel('pO_2 (PAL)','fontweight','bold');ylabel('O_3 column depth (DU)','fontweight','bold');
text(1e-3,270,'WOUDC Average = 292 DU','fontsize',14);
legend('Standard 1D','30 times lower CH3CL flux','Fixed 0.5ppt CH3CL','Fixed 0.5ppt CH3CL and No S/A','fontsize',10,'Location','best');
set(gca,'XScale','log');set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
ylim([0 350]);
hold off

subplot(212)
plot(pO2,ones(5)*2.92e8,':','Color',[0 0 0],'LineWidth',2)
hold on
plot(pO2(1:4),phi_CL,'Color',[0.93,0.69,0.13],'LineWidth',2)
plot(pO2(1:3),phi_CLNOSNOA(1:3),'Color',[0.85,0.33,0.10],'LineWidth',2)
scatter(pO2(4),phi_CLNOSNOA(4),70,'o','filled')
xlabel('pO_2 (PAL)','fontweight','bold');ylabel('CH_3CL flux (cm^{-2} s^{-1})','fontweight','bold');
%legend('Standard 1D with CH3CL','1D with fixed 0.5ppb CH3CL','1D with fixed 0.5ppb CH3CL and No S/A','Crashed run','fontsize',14,'Location','best');
set(gca,'XScale','log');set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);
%ylim([0 350]);
hold off
%%
figure
plot(pO2,mrCH3CL.*1e9,'Color',[0.4660 0.6740 0.1880],'LineWidth',2)
xlabel('pO_2 (PAL)','fontweight','bold');ylabel('CH_3Cl mixing ratio (ppb)','fontweight','bold');
set(gca,'XScale','log');%set(gca,'YScale','log');
set(gca,'LineWidth',2,'FontSize',20);