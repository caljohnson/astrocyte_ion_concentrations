function [I_current] = ncx_current(Na_out,Na_in,Ca_out,Ca_in,Vm)
%NCX_CURRENT model of the NCX (Sodium-Calcium Exchanger) current
%   relies on concentrations of sodium (Na_out, Na_in) and calcium (Ca_out, Ca_in)
%   as well as the astrocyte membrane potential (Vm) in mV

global R F T
%constants
% F = 96485; %C/mol, Faraday's constant
% R = 8.31; %J/mol K, ideal gas constant
% T = 310; %K, absolute temperature

%reversal potential in mV for the exchanger
V_NaCa = (R*T/F).*(3*log(Na_out./Na_in)-log(Ca_out./Ca_in)).*1e3; %mV

%Gall-Susa model
%hill function for Ca_in - prevents Ca_in from going negative!
nH = 1;
K_H = 1.5*1e-3; %1.5 microM
hill = Ca_in.^nH./(Ca_in.^nH + K_H.^nH);
%outputs
% I_current =  -1.*(Vm - V_NaCa);
I_current =  hill.*(Vm - V_NaCa);
% Ca_flux = -1*I_current;
% Na_flux = 3*I_current;

% %Luo-Rudy model
% KN = 87500*1e-3; %mM
% KC = 1380*1e-3; %mM
% Ksat = 0.1;
% eta = 0.35;
% 
% I_current =  (Na_out.^3)./(KN.^3 + Na_out.^3).*(Ca_out)./(KC + Ca_out)...
%     .*(Na_in.^3./(Na_out.^3).*exp(eta.*Vm.*1e-3.*F./(R*T)) - Ca_in./Ca_out.*exp((eta-1).*Vm.*1e-3.*F./(R*T)))...
%     ./(1 + Ksat.*exp((eta-1).*Vm.*1e-3.*F./(R*T)));

end

