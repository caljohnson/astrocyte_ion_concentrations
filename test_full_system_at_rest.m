%Test_Full_Syste_at_restm.m
%tracks changes in intra- and extra-cellular ion concentrations
%of extracellular ions:
%       Glu, K+, Ca2+, Na+, H+ (Glu_out, K_out, Ca_out, Na_out)
%and intracellular ions:
%       K+, Ca2+, Na+ (K_in, Ca_in, Na_in)
%through EAAT2, Kir4.1, NCX, and NKA
%also keeps track of astrocyte membrane potential Vm

addpath('./src'); close all

global g_EAAT2 g_Kir41 g_NCX rho_NKA R F T ...
    VolE VolA g_leak g_Na_leak1 g_Na_leak2 g_K_leak

%constants
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature
g_EAAT2 = 0;% not tuned yet
g_Kir41 = 5e-15;
g_NCX = 1e-9;
rho_NKA = 1e-14;
VolE = 1.41e-18;  %L, from 1.41e-3 micrometers^3, extracellular/synaptic cleft volume from Handy, Lawley, Borisyuk 2018
VolA = VolE*(2); %L, arbitrary
g_leak = 0;%1e-10;
g_Na_leak1 = 0;%1e-14; %Na in-out
g_Na_leak2 = 1e-14; %Na in-further in/other astrocytes
g_K_leak = 1e-14; %K in-further in/other astrocytes

% %resting astrocyte concentrations
Na_in = 16.6; %mM - Kirischuk et al. 2012
K_in = 140; %140mM - Kirischuk et al. 2012
% K_in = 120; %120mM - Sicca et al. 2016
% Ca_in = 73*1e-6; %mM - Kirischuk et al. 2012
Ca_in = 1.6658e-04; %mM - set so that NCX = 0 at initial concentrations

%elevated synaptic cleft (extracellular) concentrations
% % Glu_out = 0.1; %mM (elevated) - Flanagan et al 2018
% Glu_out = 0; %mM (rest) - Flanagan et al 2018
% K_out = 12; %12mM (elevated by seizures) - textbook/Flanagan et al. 2018
% % K_out = 8;
% K_out = 4; %mM (elevated by one nearby action potential) - Walz 1999
K_out = 3.5; %mM - Huguet et al. 2016
% K_out = 3; %mM (rest) - Walz 1999
Na_out = 140; %mM - Kirischuk et al. 2012
Ca_out = 2; %mM - Kirischuk et al. 2012
% H_out = (10^(-7.4))*1e3; %mM - Kirischuk et al. 2012

%concentrations from Sibille 2015
Glu_out = 0; %mM
% K_out = 2.5; %mM
% Na_out = 116; %mM
% K_in = 135; %mM
% Na_in = 12; %mM
K_neuron = 135; %mM
Na_neuron = 12; %mM
%others
% Ca_out = 2; %mM - Kirischuk et al. 2012
% Ca_in = 73e-6;%mM - Kirischuk et al. 2012 


Vm = -80; %mV, resting astrocyte membrane potential Kirischuk et al. 2012

tspan = [0,3e3];

Ca_force = @(t) 0.*t;
   
X0 = [Na_in; K_in; Ca_in; Na_out; K_out; Ca_out; Glu_out; Vm];

options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1);
[t,X] = ode15s(@(t,X) system_eqns_with_ca_signal_v2(t,X,Ca_force), tspan, X0,options);

labels = {'[Na^+]_i, mM'; '[K^+]_i, mM'; '[Ca^{2+}]_i, mM';...
    '[Na^+]_e, mM'; '[K^+]_e, mM'; '[Ca^{2+}]_e, mM'; '[Glu^-]_e, mM';...
     'V_m, mV';};
for ii=1:8
    figure(ii);
    plot(t,X(:,ii),'LineWidth',4);
    ylabel(labels(ii)); xlabel('t');
    set(gca,'FontSize',20);
    xlim([0 1e3])
end

% figure(9);
% plot(t,rho_NKA.*nka_current(X(:,1),X(:,5)),'LineWidth',2);
% hold on;
% ylabel('I_{NKA}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 1e3])
% 
% figure(10);
% plot(t,g_NCX.*ncx_current(X(:,4),X(:,1),X(:,6),X(:,3),X(:,8)),'LineWidth',2);
% hold on; ylim([-2 2].*1e-12)
% ylabel('I_{NCX}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 1e3])
% 
figure(11);
VKA = (R*T/F)*log(X(:,4)./X(:,1))*1e3; %mV
plot(t,VKA,'LineWidth',2); hold on;
ylabel('E_{Na,A}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 1e3])

% figure(12);
% plot(t,g_Kir41.*kir41_current(X(:,6),X(:,3),X(:,8)),'LineWidth',2);
% ylabel('I_{Kir}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 1e3])
% 
% figure(13);
% V_NaCa = (R*T/F).*(3*log(X(:,4)./X(:,1))-log(X(:,6)./X(:,3))).*1e3; %mV
% plot(t,V_NaCa,'LineWidth',2); hold on;
% ylabel('V_{NaCa}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 1e3])
% 
% figure(14);
% %assuming neuron has a steady concentration of [K+]i = 120 mM
% VKN = (R*T/F)*log(X(:,5)./120)*1e3; %mV
% plot(t,VKN,'LineWidth',2); hold on;
% ylabel('E_{K,N}'); xlabel('t'); set(gca,'FontSize',20);
% legend('Elevated potassium, calcium response', ...
%     'Elevated potassium, no calcium response',...
%     'Resting potassium');
