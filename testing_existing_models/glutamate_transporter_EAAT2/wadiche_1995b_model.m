% EAAT2 - Glumate Transporter model
% based on Wadiche et al. 1995b

clear; close all;

%constants
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

%-----------Wadiche et al. 1995b Data
%fit to Wadiche et al., 1995b Fig2B EAAT2 I-V data
[Vm_data, I_data] = csvimport('data/Wadiche_1995b_fig2b_10microM.csv', 'columns', [1, 2] ,'noHeader', true);
%fix positive currents to zero
I_data(I_data>0) = 0;
%ionic concentrations
%intracellular concentrations
%recorded directly from oocytes, so whatever internal concentrations are
%natural - from Gillespie 1983
Na_in = 91; %mM
K_in = 4;  %mM
H_in = 10^(-7); %mM - not sure about this
Glu_in = 0; %mM - not stated
Cl_in = 59; %mM

%extracellular concentrations
%from Wadiche 1995b:
%Recording solution (ND96) contained 96 mM NaCl, 2 mM KCl, 1 mM MgCl2, 1.8 mM CaCl2, and 5 mM HEPES (pH 7.4).
Na_out = 96; %mM
K_out = 2;  %mM
H_out = 10^(-7.4); %mM
Cl_out = 96 + 2 + 1*2 + 1.8*2; %mM
Glu_out = 0.01; %mM

%Vm for I-V curve
Vm = (-140:5:100); %mV

%functions - note RT/F has units volts
%reversal potential for [Cl-]
Vcl = (R*T/(-1*F))*log(Cl_out/Cl_in)*1e3; %mV

%model params
kglu = 1.58e-2; %mM, from wadiche_1995b_model_hillfit
% gcl = ?; %chloride conductance
%gcl is unknown x(1) below
% tau = 5; %glutamate transporter turnover rate at 0mV, 5 1/s from Wadiche 1995a Fig 7
% N = 4.1e10/(6.022e23); %number of moles of transporter expressed
% tau*N*F is unknown x(2) below
% mu = ?; %Boltzmann factor determining the voltage depen- dence of the transporter turnover rate
% mu is unknown x(3) below

% I_EAAT2 = @(x,Vm) (Glu_out./(kglu+Glu_out))*((Vm-Vcl)*x(1) - tau*N*F*exp(-x(2)*Vm));
I_EAAT2 = @(x,Vm) (Glu_out./(kglu+Glu_out))*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));
% I_EAAT2 = @(x,Vm) (Glu_out./(kglu+Glu_out))*(-x(1)*exp(-x(2)*Vm));

%my fit
my_fit = lsqcurvefit(I_EAAT2,[0.001,16,0.02],Vm_data,I_data);

figure(1);clf;
plot(Vm, I_EAAT2(my_fit,Vm), 'b-','LineWidth',3);  hold on;
xlabel('membrane potential, mV'); ylabel('current, pA');
set(gca,'FontSize',20);
xlim([-125,140]); ylim([-100,10]);
plot(Vm_data, I_data, 'bo','MarkerSize',10,'LineWidth',2,'MarkerFaceColor','b')
% legend('model', 'Wadiche et al 1995b Fig1D data')
%axes lines
line([-165,140],[0,0],'LineStyle','--', 'Color', 'k');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k');

%--- testing fit vs. more [Glu]
%30 microM
[Vm_data2, I_data2] = csvimport('data/Wadiche_1995b_fig2b_30microM.csv', 'columns', [1, 2] ,'noHeader', true);
%fix positive currents to zero
I_data2(I_data2>0) = 0;
Glu_out = 0.03; %mM
I_EAAT2 = @(x,Vm) (Glu_out./(kglu+Glu_out))*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));
plot(Vm, I_EAAT2(my_fit,Vm), 'r-','LineWidth',3); 
plot(Vm_data2, I_data2, 'r^','MarkerSize',10,'LineWidth',2)

%100 microM
[Vm_data3, I_data3] = csvimport('data/Wadiche_1995b_fig2b_100microM.csv', 'columns', [1, 2] ,'noHeader', true);
%fix positive currents to zero
I_data3(I_data3>0) = 0;
Glu_out = 0.1; %mM
I_EAAT2 = @(x,Vm) (Glu_out./(kglu+Glu_out))*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));
plot(Vm, I_EAAT2(my_fit,Vm), 'm-','LineWidth',3);  
plot(Vm_data3, I_data3, 'md','MarkerSize',10,'LineWidth',2)

%300 microM
[Vm_data4, I_data4] = csvimport('data/Wadiche_1995b_fig2b_300microM.csv', 'columns', [1, 2] ,'noHeader', true);
%fix positive currents to zero
I_data4(I_data4>0) = 0;
Glu_out = 0.3; %mM
I_EAAT2 = @(x,Vm) (Glu_out./(kglu+Glu_out))*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));
plot(Vm, I_EAAT2(my_fit,Vm), 'c-','LineWidth',3);  
plot(Vm_data4, I_data4, 'cs','MarkerSize',10,'LineWidth',2)

%1 mM
[Vm_data5, I_data5] = csvimport('data/Wadiche_1995b_fig2b_1mM.csv', 'columns', [1, 2] ,'noHeader', true);
%fix positive currents to zero
I_data5(I_data5>0) = 0;
Glu_out = 1; %mM
I_EAAT2 = @(x,Vm) (Glu_out./(kglu+Glu_out))*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));
plot(Vm, I_EAAT2(my_fit,Vm), 'k-','LineWidth',3);  
plot(Vm_data5, I_data5, 'ko','MarkerSize',10,'LineWidth',2)


%--- testing vs. 104mM and 0mM Cl current
%first fit to 104mM [Cl]o I-V data (Fig3b)
Glu_out = 1; %mM
Cl_out = 104; %mM
[Vm_data6, I_data6] = csvimport('data/Wadiche_1995b_fig3b_104mMCl.csv', 'columns', [1, 2] ,'noHeader', true);
%use fit
Vcl = (R*T/(-1*F))*log(Cl_out/Cl_in)*1e3; %mV
I_EAAT2 = @(x,Vm) (Glu_out./(kglu+Glu_out))*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));
my_fit2 = lsqcurvefit(I_EAAT2,my_fit,Vm_data6,I_data6);

%plot
Vm = (-140:5:100);
figure(2);clf;
plot(Vm, I_EAAT2(my_fit2,Vm), 'b-','LineWidth',3);  hold on;
xlabel('membrane potential, mV'); ylabel('current, pA');
set(gca,'FontSize',20);
xlim([-125,140]); ylim([-100,10]);
plot(Vm_data6, I_data6, 'bs','MarkerSize',10,'LineWidth',2)
% legend('model', 'Wadiche et al 1995b Fig1D data')
%axes lines
line([-165,140],[0,0],'LineStyle','--', 'Color', 'k');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k');

%next show against 0mM [Cl]o I-V data 
[Vm_data7, I_data7] = csvimport('data/Wadiche_1995b_fig3b_0mMCl.csv', 'columns', [1, 2] ,'noHeader', true);
Cl_out = 0; %mM
I_EAAT2_noCl = @(x,Vm) (Glu_out./(kglu+Glu_out))*(0 - x(1)*exp(-x(2)*Vm));
plot(Vm, I_EAAT2_noCl(my_fit2(2:3),Vm), 'r--','LineWidth',3);
plot(Vm_data7, I_data7, 'rs','MarkerSize',10,'LineWidth',2)


%fit just against 0mM [Cl]o
my_fit3 = lsqcurvefit(I_EAAT2_noCl,my_fit2(2:3),Vm_data7,I_data7);
figure(3); clf;
plot(Vm, I_EAAT2_noCl(my_fit3,Vm), 'r--','LineWidth',3); hold on
plot(Vm_data7, I_data7, 'rs','MarkerSize',10,'LineWidth',2)
xlabel('membrane potential, mV'); ylabel('current, pA');
set(gca,'FontSize',20);
xlim([-125,140]); ylim([-100,10]);
%axes lines
line([-165,140],[0,0],'LineStyle','--', 'Color', 'k');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k');




