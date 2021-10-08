function [I_current] = eaat2_current(Na_out,Glu_out, H_out, K_in, K_out,Vm)
%EAAT2_CURRENT model of the EAAT2 (Excitatory Amino-Acid Transporter)
% glutamate transporter current
%   relies on concentrations of external glutamate, sodium, and hydrogen 
%   (Na_out, Glu_out,H_out) and internal potassium (K_in)
%   as well as the astrocyte membrane potential (Vm) in mV

global R F T
%constants
% F = 96485; %C/mol, Faraday's constant
% R = 8.31; %J/mol K, ideal gas constant
% T = 310; %K, absolute temperature


%values from Darshan's figure
kglu = 0.2; %mM - Hill fit has degree n = 1
kNa = 350; %mM - Hill fit has degree n = 1
beta = 0.0165; %1/mV

%old values from papers
% kglu = 1.58e-2; %mM, from wadiche_1995b_model_hillfit
% x = 0.0204; %1/mV from wadiche_1995b_model_hillfit
% kNa = 46; %mM, from Zerangue & Kavanaugh 1996 Fig 3
% kH = 26e-6; %mM, from Zerangue & Kavanaugh 1996 Fig 3
% kK = 17; %mM, from Zerangue & Kavanaugh 1996 Fig 3
% kK = 12; %mM, backwards hill fit

%outputs
I_current = (Na_out./(kNa+Na_out)).*(Glu_out./(kglu+Glu_out)).*(-exp(-beta.*Vm));

%with no K_in dependence
% I_current = (Na_out.^2.3./(kNa.^2.3+Na_out.^2.3)).*(H_out./(kH+H_out)).*...
%     (Glu_out./(kglu+Glu_out)).*(exp(-x*Vm));
%with K_out dependence
% I_current = (Na_out.^2.3./(kNa.^2.3+Na_out.^2.3)).*(H_out./(kH+H_out)).*...
%     (1-0.3.*(K_out./(kK+K_out))).*(Glu_out./(kglu+Glu_out)).*(-exp(-x*Vm));

% K_flux = -1/2F*I_current;
% Glu_flux = 1/2F*I_current;
% Na_flux = 3/2F*I_current;
% H_flux = 1/2F*I_current;

end

