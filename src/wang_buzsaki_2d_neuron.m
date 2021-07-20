function [dXdt] = wang_buzsaki_2d_neuron(t,X,Iapp,ENa,EK,Vrest)
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
phi = 1;
EL = -65; %mV
% ENa = 55;%mV - i made these inputs to the code, not necessary!
% EK = -90;%mV

%functions
alpham = @(V) 0.1.*(V+35.0)./(1.0-exp(-(V+35.0)./10.0));
betam = @(V) 4.0.*exp(-(V+60.0)./18.0);
Minf = @(v) alpham(v)./(alpham(v)+betam(v));

alphah = @(V) 0.07.*exp(-(V+58.0)./20.0);
dalphahdv = @(V) (-1/20).*0.07.*exp(-(V+58.0)./20.0);
betah = @(V) 1.0./(1.0+exp(-(V+28.0)./10.0));
dbetahdv = @(V) (1/10).*exp(-(V+28.0)./10.0).*(1.0+exp(-(V+28.0)./10.0)).^(-2);
Hinf = @(v) alphah(v)./(alphah(v)+betah(v));
dHinfdv = @(v) ((alphah(v) + betah(v)).*dalphahdv(v)...
                - alphah(v).*(dalphahdv(v) + dbetahdv(v)))./ ...
                (alphah(v)+betah(v)).^2;
% tauH = @(v) 1.0./(alphah(v)+betah(v));
tauH = @(v) (1/phi).*1.0./(alphah(v)+betah(v));

alphan = @(V) 0.01.*(V+34.0)./(1.0-exp(-(V+34.0)./10.00));
dalphandv = @(V) 0.01.*((1.0-exp(-(V+34.0)./10.00) - ...
                (V+34).*(1/10).*exp(-(V+34.0)./10.00))./ ...
                (1.0-exp(-(V+34.0)./10.00)).^2 );
betan = @(V) 0.125.*exp(-(V+44.0)./80.0);
dbetandv = @(V) (-1/80).*0.125.*exp(-(V+44.0)./80.0);
Ninf = @(v) alphan(v)./(alphan(v)+betan(v));
dNinfdv = @(v) ((alphan(v) + betan(v)).*dalphandv(v)...
                - alphan(v).*(dalphandv(v) + dbetandv(v)))./ ...
                (alphan(v)+betan(v)).^2;
% tauN = @(v) 1.0./(alphan(v)+betan(v));
tauN = @(v) (1/phi).*1.0./(alphan(v)+betan(v));

% %leak current
% I_L = @(V) g_L.*(V - E_L);
% %transient sodium current
% I_Na = @(V,h) g_Na.*(m_inf(V)).^3.*h.*(V-E_Na); 
% %delayed rectifier potassium current
% I_K = @(V,n) g_K.*n.^4.*(V-E_K);

%new coordinates
% x = n - Ninf(Vrest);
% y = h - Hinf(Vrest);
%rotation angle
angle = atan(dHinfdv(Vrest)./dNinfdv(Vrest))
%projection
% z1 = cos(angle).*x + sin(angle).*y;
nnew = @(z1) Ninf(Vrest) + z1.*cos(angle);
hnew = @(z1) Hinf(Vrest) + z1.*sin(angle);
% dz1dt = cos(angle).*dndt + sin(angle).*dhdt;
dz1dt = @(V,z1) -cos(angle).*(z1.*cos(angle) + Ninf(Vrest) -Ninf(V))./tauN(V) ...
    -sin(angle).*(z1.*sin(angle) + Hinf(Vrest) - Hinf(V))./tauH(V);

%h = b - an
% hnew = @(n) 1.021 - 2.68.*n;

%state variables equations V,h,n
% dVdt = @(Ve,z1) (Iapp-gLe.*(Ve-EL)-gNa.*((Minf(Ve)).^3).*(hnew(z1)).*(Ve-ENa)...
%             -gK.*(nnew(z1).^4).*(Ve-EK))./Cm;
dVdt = @(Ve,z1) (Iapp(t)-gLe.*(Ve-EL)-gNa.*((Minf(Ve)).^3).*(hnew(z1)).*(Ve-ENa)...
            -gK.*(nnew(z1).^4).*(Ve-EK))./Cm;
% dVdt = @(Ve,n) (Iapp-gLe.*(Ve-EL)-gNa.*((Minf(Ve)).^3).*(hnew(n)).*(Ve-ENa)...
%             -gK.*(n.^4).*(Ve-EK))./Cm;
% dhdt = @(Ve,he) phi.*(Hinf(Ve)-he)./tauH(Ve);
% dndt = @(Ve,ne) phi.*(Ninf(Ve)-ne)./tauN(Ve);

%full system
dXdt = [dVdt(X(1), X(2)); dz1dt(X(1),X(2));];
% dXdt = [dVdt(X(1), X(2)); dndt(X(1),X(2));];

end

