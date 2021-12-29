function [dVdt] = EIF_neuron_ode(t,V,I)
%EIF_neuron_ode 
%   ODE for membrane potential V of an EIF neuron
%   with time-dependent input current I
%   
global V_T del_T g_L g_T E_L C

dVdt = 1/C.*(-g_L.*(V-E_L) + g_T.*del_T.*exp((V-V_T)./del_T) + I(t));

end

