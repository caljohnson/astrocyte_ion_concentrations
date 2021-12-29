%test_separate_components.m
%uses intra- and extra-cellular ion concentrations
%of extracellular ions:
%       Glu, K+, Ca2+, Na+, H+ (Glu_out, K_out, Ca_out, Na_out, H_out)
%and intracellular ions:
%       K+, Ca2+, Na+ (K_in, Ca_in, Na_in)
%to look at currents through EAAT2, Kir4.1, NCX, and NKA


addpath('./src')

%resting astrocyte concentrations
Na_in = 16.6; %mM - Kirischuk et al. 2012
K_in = 140; %mM - Kirischuk et al. 2012
Ca_in = 73*1e-6; %mM - Kirischuk et al. 2012

%elevated synaptic cleft (extracellular) concentrations
% Glu_out = 0.1; %mM (elevated) - Flanagan et al 2018
Glu_out = 0; %mM (resting)
% K_out = 12; %mM (elevated) - textbook
K_out = 3.2; %mM (resting) - model
Na_out = 140; %mM - Kirischuk et al. 2012
Ca_out = 2; %mM - Kirischuk et al. 2012
H_out = (10^(-7.4))*1e3; %mM - Kirischuk et al. 2012

Vm = -80; %mV, resting astrocyte membrane potential, Kirischuk et al. 2012

%other params
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

I_kir41 = kir41_current(K_out,K_in,Vm)
I_nka = nka_current(Na_in,K_out)
I_eaat2 = eaat2_current(Na_out,Glu_out, H_out, K_in,K_out,Vm)
V_NaCa = (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out./Ca_in))*1e3 %mV
I_ncx = ncx_current(Na_out,Na_in,Ca_out,Ca_in,Vm)
I_na_leak = na_leak_current(Na_out,Na_in,Vm)