function [dXdt] = wang_buzsaki_hippocampal_neuron_currentinject_ode(t,X,Iapp,ENa,EK)
%wang_buzsaki_hippocampal_neuron_komek_version_ode 
%   ODE for a simple Wang-Buzsaki hippocampal neuron
%   adapted from Komek et al. 2012
%   with time-dependent input current I
%   and time-dependent V_K (from astrocyte effects)

%parameters
Cm = 1; %muF/cm^2
gLe = 0.1;%mS/cm^2
gNa = 35;%mS/cm^2
gK = 9; %mS/cm^2
phi = 5;
EL = -65; %mV
% ENa = 55;%mV - i made these inputs to the code, not necessary!
% EK = -90;%mV

%functions
alpham = @(V) 0.1.*(V+35.0)./(1.0-exp(-(V+35.0)./10.0));
betam = @(V) 4.0.*exp(-(V+60.0)./18.0);
Minf = @(v) alpham(v)./(alpham(v)+betam(v));

alphah = @(V) 0.07.*exp(-(V+58.0)./20.0);
betah = @(V) 1.0./(1.0+exp(-(V+28.0)./10.0));
Hinf = @(v) alphah(v)./(alphah(v)+betah(v));
tauH = @(v) 1.0./(alphah(v)+betah(v));

alphan = @(V) 0.01.*(V+34.0)./(1.0-exp(-(V+34.0)./10.00));
betan = @(V) 0.125.*exp(-(V+44.0)./80.0);
Ninf = @(v) alphan(v)./(alphan(v)+betan(v));
tauN = @(v) 1.0./(alphan(v)+betan(v));

% %leak current
% I_L = @(V) g_L.*(V - E_L);
% %transient sodium current
% I_Na = @(V,h) g_Na.*(m_inf(V)).^3.*h.*(V-E_Na); 
% %delayed rectifier potassium current
% I_K = @(V,n) g_K.*n.^4.*(V-E_K);

%state variables equations V,h,n
dVdt = @(Ve,he,ne) (Iapp(t)-gLe.*(Ve-EL)-gNa.*((Minf(Ve)).^3).*he.*(Ve-ENa)...
            -gK.*(ne.^4).*(Ve-EK))./Cm;
dhdt = @(Ve,he,ne) phi.*(Hinf(Ve)-he)./tauH(Ve);
dndt = @(Ve,he,ne) phi.*(Ninf(Ve)-ne)./tauN(Ve);

%full system
dXdt = [dVdt(X(1), X(2), X(3)); dhdt(X(1),X(2),X(3)); dndt(X(1),X(2),X(3));];

end

