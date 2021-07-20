function [dXdt] = izhikevich_NaK_neuron_ode(t,X,I,V_K)
%izhikevich_NaK_neuron_ode 
%   ODE for a simple Izhikevich Na-K neuron
%   with time-dependent input current I
%   and time-dependent V_K (from astrocyte effects)

% V_K = interp1(V_Kt, V_K, t); % Interpolate the data set (VKt, VK) at times t


% phi = 0.04; %Hopf parameters
% phi = 1e-4;
g_Na = 20; 
g_K = 10;
g_L = 8;
V_Na = 60; %mV
V_L = -80; %mV
V1 = -10;
k1 = 15;
V2 = -25; %mV
k2 = 5; %mV
C = 1;
tau = 1;

Mss = @(V) 1./(1 + exp((V1-V)./k1));
Nss = @(V)  1./(1 + exp((V2-V)./k2));

% dVdt = (1/C).*(I - g_L.*(V-V_L) -g_Ca.*Mss(V).*(V-V_Ca) - g_K.*N.*(V-V_K));
% dNdt = (Nss(V) - N)/tau_N(V);

I_ion = @(V,N) g_Na.*Mss(V).*(V-V_Na) + g_K.*N.*(V-V_K) ...
        + g_L.*(V-V_L) ;

I_app = I(t);

dXdt =[(1/C).*(I_app - I_ion(X(1), X(2))); ...
        (Nss(X(1)) - X(2))./tau;];

end

