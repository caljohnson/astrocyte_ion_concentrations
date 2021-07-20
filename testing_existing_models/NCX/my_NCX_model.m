%my NCX model - Sodium-Calcium exchanger model
%   linear ohmic current
% based on Matsuoka & Hilgemann 1992, Hilgemann 1992 fits
% using astrocyte concentrations from Kirischuk et al. 2012


%constants
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 298.15; %K, absolute temperature (25 C in Kirischuk et al. 2012)

%cellular concentrations - Kirischuk et al. 2012
%resting concentrations
Na_in = 16.6; %mM
Ca_in = 73*1e-6; %mM
Na_out = 140; %mM
Ca_out = 2; %mM

%functions - note RT/F has units volts
%reversal potential for NCX
V_NaCa = (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out./Ca_in))*1e3; %mV

x = 0.01; %pA? 
I_NaCa = @(V) x.*(V - V_NaCa);

%Vm for I-V curve
Vm = (-150:5:100); %mV

figure(1);clf
plot(Vm,I_NaCa(Vm),'b-','LineWidth',4); hold on
line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-2,2],'LineStyle','--', 'Color', 'k','HandleVisibility','off');


%cellular concentrations - Kirischuk et al. 2012
%elevated concentrations
Na_in = 25; %mM
Ca_in = 4*1e-3; %mM
Na_out = 140; %mM
Ca_out = 2; %mM
%reversal potential for NCX
V_NaCa = (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out./Ca_in))*1e3; %mV
I_NaCa = @(V) x.*(V - V_NaCa);
plot(Vm,I_NaCa(Vm),'r-','LineWidth',4);
xlabel('V_m (mV)'); ylabel('current');set(gca,'FontSize',20);
lgd= legend('Resting','Elevated'); lgd.Location = 'northwest'; 
title('NCX, typical astrocyte concentrations');

