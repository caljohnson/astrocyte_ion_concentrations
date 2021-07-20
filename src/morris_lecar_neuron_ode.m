function [dXdt] = morris_lecar_neuron_ode(t,X,I,V_K)
%MORRIS_LECAR_NEURON_ODE 
%   ODE for a simple morris-lecar neuron
%   with time-dependent input current I
%   and time-dependent V_K (from astrocyte effects)

% V_K = interp1(V_Kt, V_K, t); % Interpolate the data set (VKt, VK) at times t


phi = 0.04; %Hopf parameters
% phi = 1e-4;
g_Ca = 4; 
g_K = 8;
g_L = 2;
V_Ca = 120; %mV
V_L = -60; %mV
V1 = -1.2;
V2 = 18;
V3 = 2; %mV
V4 = 30; %mV
C = 20;
% tau = 0.8;

Mss = @(V) (1 + tanh((V-V1)/V2) )/2;
Nss = @(V) (1 + tanh((V-V3)/V4) )/2;
r = @(V) cosh((V-V3)/(2*V4)) * phi;

% dVdt = (1/C).*(I - g_L.*(V-V_L) -g_Ca.*Mss(V).*(V-V_Ca) - g_K.*N.*(V-V_K));
% dNdt = (Nss(V) - N)/tau_N(V);

I_ion = @(V,N) g_Ca.*Mss(V).*(V-V_Ca) + g_K.*N.*(V-V_K) ...
        + g_L.*(V-V_L) ;

I_app = I(t);

dXdt =[(1/C).*(I_app - I_ion(X(1), X(2))); ...
        (Nss(X(1)) - X(2)).*r(X(1));];

end

