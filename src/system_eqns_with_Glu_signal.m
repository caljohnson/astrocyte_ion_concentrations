function [dX_dt] = system_eqns_with_Glu_signal(t,X,Glu_forcing,Ca_forcing)
%SYSTEM_EQNS creates system of differential equations for the
%            astrocyte ion concentration tracking model
%            with a prescribed Glu and Ca FORCING TERMS
%       dX/dt = F(X)
% where X = (Na_in, K_in, Ca_in, Na_out, K_out, Ca_out, Glu_out, H_out, Vm) 
%tracks changes in intra- and extra-cellular ion concentrations
%of extracellular ions:
%       Glu, K+, Ca2+, Na+, H+ (Glu_out, K_out, Ca_out, Na_out, H_out)
%and intracellular ions:
%       K+, Ca2+, Na+ (K_in, Ca_in, Na_in)
%through EAAT2, Kir4.1, NCX, and NKA
%also keeps track of astrocyte membrane potential Vm

global F g_EAAT2 g_Kir41 g_NCX rho_NKA VolE VolA ...
    g_leak g_Na_leak1 g_Na_leak2 g_K_leak g_Ca_leak

% g_EAAT2 = 8.381e4; %pA, from fit to Wadiche et al 1995b data
% g_EAAT2 = g_EAAT2*1e-12;%pA to C/s
% g_Kir41 = 1e5*0.05;%nA/mV, from average fit, set to Seifart et al 2009 scale
% g_Kir41 = 10*g_Kir41*1e-9;%nA/mV to C/mV/s
% g_NCX = 234;%pS, pA/V from Gall & Susa 1999 (234-1000 pS)
% g_NCX = g_NCX*1e-12*1e3; %pA/V to A/V to A/mV
% rho_NKA = 20; %muA, from Huguet et al 2015
% rho_NKA = rho_NKA*1e-6; %muA to A
% g_EAAT2 = 0;%1e-15;
% g_Kir41 = 5e-15;
% g_NCX = 1e-9;
% rho_NKA = 1e-13;
% g_Na_leak = 0;%1e-19;

%other parameters
C = 15e-12; %15pF from Sibille et al 2015
C = C*1e3; %astrocyte capacitance times 1e3 for unit conversion from V/C to mV/C
% F = 96485; %Faraday constant (96485 Coulombs/ mol)
% VolE = 1.41e-18; %L, from 1.41e-3 micrometers^3, extracellular/synaptic cleft volume from Handy, Lawley, Borisyuk 2018
% VolE = 4.16e-13; %L, from 416 micrometers^3, extracellular/synaptic cleft volume from Huguet et al. 2016
% VolA = 0.5*VolE; %astrocyte volume
% VolA = 2e-12; %L, from 2000 micrometers^3, astrocyte volume from Huguet et al. 2016
Vleak = -80; %resting membrane potential roughly -80mV, Kirischuk et al. 2012
% g_leak = 0;%1e-14;
H_out = (10^(-7.4))*1e3; %mM - Kirischuk et al. 2012
K_in_rest = 140; %mM
Na_in_rest = 16.6; %mM
Ca_in_rest = 73*1e-6; %mM - Kirischuk et al. 2012

%differential equations for the system

dVm_dt = @(Na_in, K_in, Ca_in, Na_out, K_out, Ca_out, Glu_out, Vm) ...
        -(1/C).*(rho_NKA.*nka_current(Na_in,K_out) ...
        + g_NCX.*ncx_current(Na_out,Na_in,Ca_out,Ca_in,Vm)...
        + g_Kir41.*kir41_current(K_out,K_in,Vm)...
        + g_EAAT2.*eaat2_current(Na_out,Glu_out, H_out, K_in,K_out,Vm)...
        + g_leak.*(Vm - Vleak) ...
        + g_Na_leak1.*na_leak_current(Na_out,Na_in,Vm))...
        +g_Na_leak2.*F.*VolA.*(Na_in - K_in_rest) ...
        +g_K_leak.*F.*VolA.*(K_in - K_in_rest) ...
        +g_Ca_leak.*F.*VolA.*(Ca_in - Ca_in_rest);
    
dNa_in_dt = @(Na_in, K_in, Ca_in, Na_out, K_out, Ca_out, Glu_out, Vm) ...
    (1/(VolA*F)).*(-3.*g_NCX.*ncx_current(Na_out,Na_in,Ca_out,Ca_in,Vm)...
            -3.*rho_NKA.*nka_current(Na_in,K_out) ...
            -3/2.*g_EAAT2.*eaat2_current(Na_out,Glu_out, H_out, K_in,K_out,Vm)...
            -g_Na_leak1.*na_leak_current(Na_out,Na_in,Vm))...
            -g_Na_leak2.*(Na_in - Na_in_rest);
        
dNa_out_dt = @(Na_in, K_in, Ca_in, Na_out, K_out, Ca_out, Glu_out, Vm) ...
    (1/(VolE*F)).*(+3.*g_NCX.*ncx_current(Na_out,Na_in,Ca_out,Ca_in,Vm)...
            +3.*rho_NKA.*nka_current(Na_in,K_out) ...
            +3/2.*g_EAAT2.*eaat2_current(Na_out,Glu_out, H_out, K_in,K_out,Vm)...
            +g_Na_leak1.*na_leak_current(Na_out,Na_in,Vm));
     %-(VolA/VolE).*dNa_in_dt(Na_in, K_in, Ca_in, Na_out, K_out, Ca_out, Glu_out, Vm);
        

dK_in_dt = @(Na_in, K_in, Ca_in, Na_out, K_out, Ca_out, Glu_out, Vm) ...
    (1/(VolA*F)).*(-1.*g_Kir41.*kir41_current(K_out,K_in,Vm)...
            +1.*rho_NKA.*nka_current(Na_in,K_out) ...
            +1/2.*g_EAAT2.*eaat2_current(Na_out,Glu_out, H_out, K_in,K_out,Vm))...
            -g_K_leak.*(K_in - K_in_rest);

dK_out_dt = @(Na_in, K_in, Ca_in, Na_out, K_out, Ca_out, Glu_out, Vm) ...
...%     -(VolA/VolE).*dK_in_dt(Na_in, K_in, Ca_in, Na_out, K_out, Ca_out, Glu_out, Vm);
 (1/(VolE*F)).*(+1.*g_Kir41.*kir41_current(K_out,K_in,Vm)...
            -1.*rho_NKA.*nka_current(Na_in,K_out) ...
            -1/2.*g_EAAT2.*eaat2_current(Na_out,Glu_out, H_out, K_in,K_out,Vm));


dCa_in_dt = @(Na_in, K_in, Ca_in, Na_out, K_out, Ca_out, Glu_out, Vm, Ca_forcing) ...
    (1/(VolA*F)).*(1.*g_NCX.*ncx_current(Na_out,Na_in,Ca_out,Ca_in,Vm))...
    +Ca_forcing ...
    -g_Ca_leak.*(Ca_in - Ca_in_rest); %leak back into ER

dCa_out_dt = @(Na_in, K_in, Ca_in, Na_out, K_out, Ca_out, Glu_out,  Vm) ...
    (1/(VolE*F)).*(-1.*g_NCX.*ncx_current(Na_out,Na_in,Ca_out,Ca_in,Vm));
%     -(VolA/VolE).*dCa_in_dt(Na_in, K_in, Ca_in, Na_out, K_out, Ca_out, Glu_out, Vm);
    

dGlu_out_dt = @(Na_in, K_in, Ca_in, Na_out, K_out, Ca_out, Glu_out, Vm,Glu_forcing) ...
    (1/(VolE*F)).*(1/2.*g_EAAT2.*eaat2_current(Na_out,Glu_out, H_out, K_in,K_out,Vm))...
    + Glu_forcing;

%system form
%X = (Na_in, K_in, Ca_in, Na_out, K_out, Ca_out, Glu_out,Vm)
dX_dt = ...
    [dNa_in_dt(X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8)); ...
    dK_in_dt(X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8)); ...
    dCa_in_dt(X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8), Ca_forcing(t)); ...
    dNa_out_dt(X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8)); ...
    dK_out_dt(X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8)); ...
    dCa_out_dt(X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8)); ...
    dGlu_out_dt(X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8),Glu_forcing(t)); ...
    dVm_dt(X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8));];


%     dH_out_dt(X(1), X(2), X(3), X(4), X(5), X(6), X(7), X(8)); ...

end


