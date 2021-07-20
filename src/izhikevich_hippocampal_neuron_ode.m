function [dXdt] = izhikevich_hippocampal_neuron_ode(t,X,I,V_K)
%izhikevich_hippocampal_neuron_ode 
%   ODE for a simple Izhikevich hippocampal neuron
%   with time-dependent input current I
%   and time-dependent V_K (from astrocyte effects)



I_app = I(t);

dXdt = [(1/50).*(0.5.*(X(1)+60).*(X(1)+45)-X(2)+I_app);...
        0.02.*(0.5.*(X(1)+60)-X(2));];

end

