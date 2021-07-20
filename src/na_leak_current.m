function [I_current] = na_leak_current(Na_out,Na_in,Vm)
%Na_LEAK_CURRENT model of the sodium leak current
%   relies on concentrations of sodium (Na_out, Na_in)
%   as well as the astrocyte membrane potential (Vm) in mV

global R F T
%constants
% F = 96485; %C/mol, Faraday's constant
% R = 8.31; %J/mol K, ideal gas constant
% T = 310; %K, absolute temperature

% phi = -(Vm.*1e-3.*F)./(R*T); 
E_Na = (R*T/F).*log(Na_out./Na_in).*1e3;

%outputs
I_current = (Vm - E_Na);
% I_current = phi.*(Na_out.*exp(-phi) - Na_in)./(exp(-phi)-1);

end

