% EAAT2 - Glumate Transporter model
% based on I-V data

clear; close all; clc;
%Kir
%constants
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

%-----------Levy et al., 1998 Data
%ionic concentrations - data from Levy et al. 1998
%intracellular concentrations
Na_in = 12; %mM
K_in = 140;  %mM
H_in = 10^(-7); %mM
Glu_in = 0.2; %mM
%extracellular concentrations
Na_out = 142; %mM
K_out = 2.5;  %mM
H_out = 10^(-7.4); %mM


%--------------- Fitting model to Extracellular Glutamate concentration ----
[Glu_out_data, I_Glu_data] = csvimport('data/levy_etal_1998_fig2b.csv', 'columns', [1, 2] ,'noHeader', true);
Glu_out_data = Glu_out_data*1e-3; %convert from microM to mM
Glu_out = logspace(log10(min(Glu_out_data)),log10(max(Glu_out_data)),20);
Vm = -110;
%new V_reversal, glutamate-dependent
V_rev = (R*T/(2*F)).*log((Na_out/Na_in)^3*(K_in/K_out)*(H_out/H_in).*(Glu_out./Glu_in))*1e3; %mV

y = [.05, -0.0134, 61.18];
%remake my model function with new Glu_out since it depends on it
I_EAAT2 = @(x,Glu) y(1).*(Vm-y(3)).*exp(y(2).*(Vm-y(3))).*(Glu.^x(1)./(x(2).^x(1)+Glu.^x(1)));
I_norm = @(x,Glu) -I_EAAT2(x,Glu)./max(abs(I_EAAT2(x,Glu)));

%fit my model
my_fit = lsqcurvefit(I_norm,[2,1.5e-2],Glu_out_data,I_Glu_data);

%plot
figure(3);clf;
semilogx(Glu_out, -I_EAAT2(my_fit,Glu_out)./max(abs(I_EAAT2(my_fit,Glu_out))), 'r-','LineWidth',3);
hold on;
xlabel('extracellular Glutamate concentration, mM'); ylabel('current, normalized');
set(gca,'FontSize',20);
semilogx(Glu_out_data, I_Glu_data, 'm^','MarkerSize',10,'LineWidth',2)
legend('my model', 'Levy et al., 1998 data')

%------ Testing vs. I-V data
Glu_out = 0.1; %mM
%Vm for I-V curve
Vm = (-125:5:140); %mV

%functions - note RT/F has units volts
V_rev = (R*T/(2*F))*log((Na_out/Na_in)^3*(K_in/K_out)*(H_out/H_in)*(Glu_out/Glu_in))*1e3; %mV

%my model
kglu = my_fit(1);
nglu = my_fit(2);
I_EAAT2 = @(Vm) y(1).*(Vm-y(3)).*exp(y(2).*(Vm-y(3))).*(Glu_out.^nglu./(kglu.^nglu+Glu_out.^nglu));

%compare with Levy et al., 1998 Fig1 (B) Glutate-transporter GLT-1 (EAAT2 analog) I-V data
[Vm_data, I_data] = csvimport('data/levy_etal_1998_fig1b.csv', 'columns', [1, 2] ,'noHeader', true);
%add (V_rev, 0) to dataset
Vm_data(end+1) = V_rev;
I_data(end+1) = 0;

figure(1);clf;
plot(Vm, I_EAAT2(Vm), 'r-','LineWidth',3); hold on
xlabel('membrane potential, mV'); ylabel('current, pA');
set(gca,'FontSize',20);
xlim([-125,140]); ylim([-100,0]);
plot(Vm_data, I_data, 'm^','MarkerSize',10,'LineWidth',2)
legend('my model', 'Levy et al., 1998 data')
%axes lines
line([-165,140],[0,0],'LineStyle','--', 'Color', 'k');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k');


%------------Bergles and Jahr, 1997 Data
%---- Fitting models to Extracellular Glutamate w/ Bergles & Jahr data
[Glu_out_data, I_Glu_data] = csvimport('data/bergles_jahr_1997_fig8c_peak.csv', 'columns', [1, 2] ,'noHeader', true);
Glu_out = logspace(log10(min(Glu_out_data)),log10(max(Glu_out_data)),20);
Vm = -100;
%new V_reversal, glutamate-dependent
V_rev = (R*T/(2*F)).*log((Na_out/Na_in)^3*(K_in/K_out)*(H_out/H_in).*(Glu_out./Glu_in))*1e3; %mV

%remake my model function with new Glu_out since it depends on it
% I_EAAT2 = @(x,Glu) y(1).*(Vm-y(3)).*exp(y(2).*(Vm-y(3))).*(Glu.^x(1)./(x(2).^x(1)+Glu.^x(1)));
I_EAAT2 = @(x,Glu) (Glu.^x(1)./(x(2).^x(1)+Glu.^x(1)));
I_norm = @(x,Glu) -I_EAAT2(x,Glu)./max(abs(I_EAAT2(x,Glu)));

%fit my model
my_fit2 = lsqcurvefit(I_norm,[1.1,1.5e-2],Glu_out_data,I_Glu_data);

%plot
figure(4);clf;
semilogx(Glu_out, -I_EAAT2(my_fit2,Glu_out)./max(abs(I_EAAT2(my_fit2,Glu_out))), 'r-','LineWidth',3);
hold on
xlabel('extracellular Glutamate concentration, mM'); ylabel('current, normalized');
set(gca,'FontSize',20);
semilogx(Glu_out_data, I_Glu_data, 'm^','MarkerSize',10,'LineWidth',2)
legend('my model', ...
   'Bergles and Jahr, 1997 data')

%fit to Bergles and Jahr, 1997 Fig7(B) Glutate-transporter EAAT2 I-V data
[Vm_data2, I_data2] = csvimport('data/bergles_jahr_1997_fig7b.csv', 'columns', [1, 2] ,'noHeader', true);
%scale current data to other paper's data
I_data2 = I_data2.*max(abs(I_data));

%Vm for I-V curve
Vm = (-125:5:140); %mV

%intracellular concentrations
Na_in = 0.01; %mM %total guess based on other data, not reported in paper
K_in = 130;  %mM
H_in = 10^(-7.2); %mM
Glu_in = 1; %mM %total guess based on other data, not reported in paper
%extracellular concentrations
Na_out = 135; %mm
K_out = 5.4;  %mM %typo in Flanagan's SI, they say they used 4.5mM but Bergles & Jahr report 5.4mM
H_out = 10^(-7.2); %mM
Glu_out = 10; %mM %also a guess based on the paper mentioning 10mM glutamate used to activate the channel

%back to single V_reversal
V_rev = (R*T/(2*F))*log((Na_out/Na_in)^3*(K_in/K_out)*(H_out/H_in)*(Glu_out/Glu_in))*1e3; %mV

%my model
kglu = my_fit(1);
nglu = my_fit(2);
I_EAAT2 = @(Vm) y(1).*(Vm-y(3)).*exp(y(2).*(Vm-y(3))).*(Glu_out.^nglu./(kglu.^nglu+Glu_out.^nglu));

%add (V_rev, 0) to dataset
Vm_data2(end+1) = V_rev;
I_data2(end+1) = 0;


figure(2);clf;
plot(Vm, I_EAAT2(Vm), 'r-','LineWidth',3);  hold on;
xlabel('membrane potential, mV'); ylabel('current, pA');
set(gca,'FontSize',20);
xlim([-125,140]); ylim([-100,0]);
plot(Vm_data2, I_data2, 'm^','MarkerSize',10,'LineWidth',2)
legend('my model', 'Bergles and Jahr, 1997 data')
%axes lines
line([-165,140],[0,0],'LineStyle','--', 'Color', 'k');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k');
