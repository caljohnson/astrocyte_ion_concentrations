%astrocyte_modulatory_role_elevatedK_tuning.m
%First, we balanced the component contributions to the full system so 
%that the system is at rest at relevant physiological ionic concentrations. 
%In this process, we tuned the contributions of NCX, NKA, and Kir4.1, 
%as well as individual leak currents. 
%To aid in this process, we also fit the elevated-potassium clearance of 
%[16] by tuning the relative strengths of Kir4.1 and NKA in the model,
%so that elevated external potassium was returned to baseline by the 
%astrocyte on the order of tens of seconds (roughly 20 s in our model).

%tracks changes in intra- and extra-cellular ion concentrations
%of extracellular ions:
%       Glu, K+, Ca2+, Na+, H+ (Glu_out, K_out, Ca_out, Na_out)
%and intracellular ions:
%       K+, Ca2+, Na+ (K_in, Ca_in, Na_in)
%through EAAT2, Kir4.1, NCX, and NKA
%also keeps track of astrocyte membrane potential Vm

addpath('./src'); close all

global g_EAAT2 g_Kir41 g_NCX rho_NKA R F T ...
    VolE VolA g_leak g_Na_leak1 g_Na_leak2 g_K_leak g_Ca_leak

%constants
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature
% g_EAAT2 = 5e-14;
g_EAAT2 = 1.2e-13; %new - 10/19
% g_Kir41 = 1e-14;
g_Kir41 = 1e-13; %new - 10/19
% g_NCX = 1e-9;
g_NCX = 5e-15; %new - 10/19
% rho_NKA = 2e-14;
rho_NKA = 24e-14; %new - 10/19
VolE = 1.41e-18;  %L, from 1.41e-3 micrometers^3, extracellular/synaptic cleft volume from Handy, Lawley, Borisyuk 2018
VolA = VolE*(1); %L, arbitrary
g_leak = 0;%1e-14;
% g_Na_leak1 = 2.1848e-16; %Na in<->out
g_Na_leak1 = 2.6217e-15; %new 10/19
g_Na_leak2 = 0;%1e-14; %Na in->further in/other astrocytes
g_K_leak = 3.68e-2; %K in->further in/other astrocytes
g_Ca_leak = 1e2; %Ca in -> ER

%resting astrocyte concentrations
Na_in = 16.6; %mM - Kirischuk et al. 2012
K_in = 140; %140mM - Kirischuk et al. 2012
Ca_in = 73*1e-6; %mM - Kirischuk et al. 2012

%elevated synaptic cleft (extracellular) concentrations
Glu_out = 0; %mM (rest)
K_out = 9; %12mM (elevated by seizures) - textbook/Flanagan et al. 2018
Na_out = 140; %mM - Kirischuk et al. 2012

Ca_out = 0.87; %mM - rest state in this model
H_out = (10^(-7.4))*1e3; %mM - Kirischuk et al. 2012

Vm = -80; %mV, resting astrocyte membrane potential Kirischuk et al. 2012

Na_outs = {};
K_outs = {};
ts = {};
tspan = [0,1e3];
tmaxplot = 5e2;

for jj=1:2
if jj==1
    Ca_force = @(t) 0.*t;
else
    Ca_force = @(t) 0.*t;
    K_out = X(end,5); %3.2;
end
    
X0 = [Na_in; K_in; Ca_in; Na_out; K_out; Ca_out; Glu_out; Vm];

options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1);
[t,X] = ode15s(@(t,X) system_eqns_with_ca_signal_v2(t,X,Ca_force), tspan, X0,options);
dt = 1e-1;
t0 = tspan(1):dt:tspan(2);
X = interp1(t,X,t0);
t = t0;

labels = {'[Na^+]_i, mM'; '[K^+]_i, mM'; '[Ca^{2+}]_i, mM';...
    '[Na^+]_e, mM'; '[K^+]_e, mM'; '[Ca^{2+}]_e, mM'; '[Glu^-]_e, mM';...
     'V_m, mV';};
for ii=1:8
    figure(ii);
    
    if jj==1
        plot(t,X(:,ii),'-','LineWidth',4);
    else
        plot(t,X(:,ii),'--','LineWidth',4);
    end
    hold on;
    
    ylabel(labels(ii)); xlabel('t (sec)');
    set(gca,'FontSize',20);
    legend('Elevated potassium',...
    'Resting potassium');
    xlim([0 tmaxplot])
end

figure(9);
if jj==1
        plot(t,rho_NKA.*nka_current(X(:,1),X(:,5)),'LineWidth',4);
    elseif jj==2
        plot(t,rho_NKA.*nka_current(X(:,1),X(:,5)),'--','LineWidth',4);
    else
        plot(t,rho_NKA.*nka_current(X(:,1),X(:,5)),'-.','LineWidth',4);
end
hold on;
ylabel('I_{NKA}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 tmaxplot])
legend('Elevated potassium, calcium response', ...
    'Elevated potassium, no calcium response',...
    'Resting potassium');

figure(10);
if jj==1
        plot(t,g_NCX.*ncx_current(X(:,4),X(:,1),X(:,6),X(:,3),X(:,8)),'LineWidth',4);
    elseif jj==2
        plot(t,g_NCX.*ncx_current(X(:,4),X(:,1),X(:,6),X(:,3),X(:,8)),'--','LineWidth',4);
    else
        plot(t,g_NCX.*ncx_current(X(:,4),X(:,1),X(:,6),X(:,3),X(:,8)),'-.','LineWidth',4);
end
hold on; ylim([-2 2].*1e-12)
ylabel('I_{NCX}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 tmaxplot])
legend('Elevated potassium, calcium response', ...
    'Elevated potassium, no calcium response',...
    'Resting potassium');

% figure(11);
% VKA = (R*T/F)*log(X(:,5)./X(:,2))*1e3; %mV
% plot(t,VKA,'LineWidth',2); hold on;
% ylabel('E_{K,A}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 tmaxplot])

figure(12);
    if jj==1
        plot(t,g_Kir41.*kir41_current(X(:,5),X(:,2),X(:,8)),'LineWidth',4);
    elseif jj==2
        plot(t,g_Kir41.*kir41_current(X(:,5),X(:,2),X(:,8)),'--','LineWidth',4);
    else
        plot(t,g_Kir41.*kir41_current(X(:,5),X(:,2),X(:,8)),'-.','LineWidth',4);
    end
hold on
ylabel('I_{Kir}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 tmaxplot])
legend('Elevated potassium, calcium response', ...
    'Elevated potassium, no calcium response',...
    'Resting potassium');

figure(13);
if jj==1
        plot(t,g_Na_leak1.*na_leak_current(X(:,4),X(:,1), X(:,8)),'LineWidth',4);
    elseif jj==2
        plot(t,g_Na_leak1.*na_leak_current(X(:,4),X(:,1), X(:,8)),'--','LineWidth',4);
    else
        plot(t,g_Na_leak1.*na_leak_current(X(:,4),X(:,1), X(:,8)),'-.','LineWidth',4);
end
hold on
ylabel('I_{NaLeak}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 tmaxplot])
legend('Elevated potassium, calcium response', ...
    'Elevated potassium, no calcium response',...
    'Resting potassium');


Na_outs{jj} = X(:,4);
K_outs{jj} = X(:,5);
ts{jj} = t;

save('Na_K_outs_elevatedK_noCa.mat','ts','Na_outs', 'K_outs');
end
figure(5);
