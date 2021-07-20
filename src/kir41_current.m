function [I_current] = kir41_current(K_out,K_in,Vm)
%KIR41_CURRENT model of the kir4.1 (Inward-Rectifying Potassium channel) current
%   relies on concentrations of potassium (K_out, K_in)
%   as well as the astrocyte membrane potential (Vm) in mV
global R F T
%constants
% F = 96485; %C/mol, Faraday's constant
% R = 8.31; %J/mol K, ideal gas constant
% T = 310; %K, absolute temperature

%reversal potential in mV for the potassium channel
VKA = (R*T/F)*log(K_out./K_in)*1e3; %mV

%using Chai fit values
% x = [0.0482; -0.1059; 269.7807;];
% x = [0.0482; 0; 269.7807;]; %with no shift at all
% x = [0.0482; +10; 269.7807;]; %with big shift
%using Sicca Fit values
% x = [0.0902;   13.7412;   85.1605]
% x = [0.0459;  17.8060;   85.1607;];
%using Sicca fit + adjusted so that I = 0 at K_out = 3.5, K_in = 140 mM
 x = [0.0459;  18.4912;   85.1607;];

%outputs
I_current =  (Vm-VKA-x(2))./(1+exp((Vm-VKA-x(2))./x(3))); %units nA
% K_flux = -1*I_current;

end

