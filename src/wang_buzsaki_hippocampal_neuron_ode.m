function [dXdt] = wang_buzsaki_hippocampal_neuron_ode(t,X,I_app,E_Na,E_K)
%wang_buzsaki_hippocampal_neuron_ode 
%   ODE for a simple Wang-Buzsaki hippocampal neuron
%   with time-dependent input current I
%   and time-dependent V_K (from astrocyte effects)

%parameters
phi = 5;
Cm = 1; %muF/cm^2

%other params
g_L = 0.1;%mS/cm^2
g_Na = 35;%mS/cm^2
g_K = 9; %mS/cm^2
E_L = -65; %mV
% E_Na = 55;%mV - i made these inputs to the code, not necessary!
% E_K = -90;%mV

%functions
alpha_m = @(V) -0.1.*(V+35)./(exp(-0.1.*(V+35))-1);
beta_m = @(V) 4.*exp(-(V+60)./18);
m_inf = @(V) alpha_m(V)./(alpha_m(V) + beta_m(V));

alpha_h = @(V) 0.07.*exp(-(V+58)./20);
beta_h = @(V) 1./(exp(-0.1.*(V+28))+1);
h_inf = @(V) alpha_h(V)./(alpha_h(V) + beta_h(V));
tau_h = @(V) 1.0./(alpha_h(V)+beta_h(V));

alpha_n = @(V) -0.01.*(V+34)./(exp(-0.1.*(V+34))-1);
beta_n = @(V) 0.125.*(exp(-(V+44)./80));
n_inf = @(V) alpha_n(V)./(alpha_n(V) + beta_n(V));
tau_n = @(V) 1.0./(alpha_n(V)+beta_n(V));

%leak current
I_L = @(V) g_L.*(V - E_L);
%transient sodium current
I_Na = @(V,h) g_Na.*(m_inf(V)).^3.*h.*(V-E_Na); 
%delayed rectifier potassium current
I_K = @(V,n) g_K.*n.^4.*(V-E_K);

%applied current at time t
% I_app = 0;% I(t);

%state variables V,h,n
dVdt = @(V,h,n) (1/Cm).*(I_app -I_Na(V,h) -I_K(V,n) -I_L(V) );
dhdt = @(V,h,n) phi.*(h_inf(V)-h)./tau_h(V); %these are the version from Komek et al. 2012
dndt = @(V,h,n) phi.*(n_inf(V)-n)./tau_n(V); %these are the version from Komek et al. 2012
% dhdt = @(V,h,n) phi.*(alpha_h(V).*(1-h) - beta_h(V).*h); %original equations from Wang-Buzsaki
% dndt = @(V,h,n) phi.*(alpha_n(V).*(1-n) - beta_n(V).*n); %original equations from Wang-Buzsaki

%full system
dXdt = [dVdt(X(1), X(2), X(3)); dhdt(X(1),X(2),X(3)); dndt(X(1),X(2),X(3));];

end

