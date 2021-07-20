%Test_Full_System.m
%tracks changes in intra- and extra-cellular ion concentrations
%of extracellular ions:
%       K+, Ca2+, Na+, H+ (K_out, Ca_out, Na_out)
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
g_EAAT2 = 3e-14;
g_Kir41 = 1e-14;
g_NCX = 1e-9;
rho_NKA = 1e-13;
VolE = 1.41e-18;  %L, from 1.41e-3 micrometers^3, extracellular/synaptic cleft volume from Handy, Lawley, Borisyuk 2018
% VolA = VolE*(1/2); %L, arbitrary - this one gives the right internal Na/Ca shapes with Glu=2mM
VolA = VolE*(1); %L, arbitrary - this is what they use in Oschmann et al. 2017
g_leak = 0;%1e-14;
g_Na_leak1 = 2.65e-16; %Na in<->out
g_Na_leak2 = 0;%1e-14; %Na in->further in/other astrocytes
g_K_leak = 3.68e-2; %K in->further in/other astrocytes
g_Ca_leak = 1e2; %Ca in -> ER

%resting astrocyte concentrations
Na_in = 16.6; %mM - Kirischuk et al. 2012
K_in = 140; %140mM - Kirischuk et al. 2012
% K_in = 120; %120mM - Sicca et al. 2016
Ca_in = 73*1e-6; %mM - Kirischuk et al. 2012
% Ca_in = 1.6658e-04; %mM - set so that NCX = 0 at initial concentrations

%elevated synaptic cleft (extracellular) concentrations
% Glu_out = 0.1; %mM (elevated) - Flanagan et al 2018
% Glu_out = 2; %mM (elevated) - Kirischuk et al. 2007 -  this one gives the right internal Na/Ca shapes w/ VolA=1/2VolE
% Glu_out = 10; %mM - bath amount
% Glu_out = 0; %mM (rest) - Flanagan et al 2018
% K_out = 12; %12mM (elevated by seizures) - textbook/Flanagan et al. 2018
% K_out = 5; %Kirischuk et al. 2012
% K_out = 3.5; %Huguet et al. 2016
% K_out = 4; %mM (elevated by one nearby action potential) - Walz 1999
% K_out = 3; %mM (rest) - Walz 1999
K_out = 3.2; %rest in this system
Na_out = 140; %mM - Kirischuk et al. 2012
% Ca_out = 2; %mM - Kirischuk et al. 2012
Ca_out = 0.9; %mM - rest state in this model
H_out = (10^(-7.4))*1e3; %mM - Kirischuk et al. 2012

Vm = -80; %mV, resting astrocyte membrane potential Kirischuk et al. 2012

% X0 = [Na_in; K_in; Ca_in; Na_out; K_out; Ca_out; Glu_out; Vm];

Na_outs = {};
K_outs = {};
ts = {};
tspan = [0,5e2];
tmaxplot = 5e2;

Glu_bath = @(t) 0*t + 1*(t<20); %1mM bath for 20 seconds

for jj=1:2
if jj==1
    Ca_force = @(t) 0.*t + 1e-1.*(t-0).*(t>=0 & t<5) + ...
        (5e-1 - 1e-1.*(t-5)).*(t>=5 & t<10) ; %one peak
% elseif jj==2
%      Ca_force = @(t) 0.*t + 1e-13.*(t-0).*(t>=0 & t<5) + ...
%         (5e-13 - 1e-13.*(t-5)).*(t>=5 & t<10) +...
%         0.5e-13.*(t-25).*(t>=25 & t<30) + ...
%         (2.5e-13 - 0.5e-13.*(t-30)).*(t>=30 & t<35); %two peak
% elseif jj==3
%    Ca_force = @(t) 0.*t + 1e-13.*(t-0).*(t>=0 & t<5) + ...
%         (5e-13 - 1e-13.*(t-5)).*(t>=5 & t<10) +...
%         0.3e-13.*(t-15).*(t>=15 & t<20) + ...
%         (1.5e-13 - 0.3e-13.*(t-20)).*(t>=20 & t<25); %sustained peak
elseif jj==2
    Ca_force = @(t) 0.*t;
end
   
X0 = [Na_in; K_in; Ca_in; Na_out; K_out; Ca_out; Vm];

options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1);
[t,X] = ode15s(@(t,X) system_eqns_with_Glu_bath(t,X,Glu_bath,Ca_force), tspan, X0,options);
dt = 1e-1;
t0 = tspan(1):dt:tspan(2);
X = interp1(t,X,t0);
t = t0;

labels = {'[Na^+]_i, mM'; '[K^+]_i, mM'; '[Ca^{2+}]_i, mM';...
    '[Na^+]_e, mM'; '[K^+]_e, mM'; '[Ca^{2+}]_e, mM';...
     'V_m, mV';};
for ii=1:7
    figure(ii);
    
    if jj==1
        plot(t,X(:,ii),'-','LineWidth',4);
    elseif jj==2
        plot(t,X(:,ii),'--','LineWidth',4);
    else
        plot(t,X(:,ii),'g-.','LineWidth',4);
    end
    hold on;
    
    ylabel(labels(ii)); xlabel('t');
    set(gca,'FontSize',20);
    legend('Glu bath, calcium transient', ...
    'Glu bath, no calcium transient');
    xlim([0 tmaxplot])
    if ii==7
        ylim([-100,-60]);
    end
end

% figure(9);
% Ca_forcing = Ca_force(t);
% plot(t,Ca_forcing,'LineWidth',2); hold on;
%     ylabel('prescribed [Ca^{2+}]_i forcing term'); xlabel('t');
%     set(gca,'FontSize',20);

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
    legend('Glu bath, calcium transient', ...
    'Glu bath, no calcium transient');

figure(10);
if jj==1
        plot(t,g_NCX.*ncx_current(X(:,4),X(:,1),X(:,6),X(:,3),X(:,7)),'LineWidth',4);
    elseif jj==2
        plot(t,g_NCX.*ncx_current(X(:,4),X(:,1),X(:,6),X(:,3),X(:,7)),'--','LineWidth',4);
    else
        plot(t,g_NCX.*ncx_current(X(:,4),X(:,1),X(:,6),X(:,3),X(:,7)),'-.','LineWidth',4);
end
hold on; ylim([-2 2].*1e-12)
ylabel('I_{NCX}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 tmaxplot])
    legend('Glu bath, calcium transient', ...
    'Glu bath, no calcium transient');

figure(11);
plot(t,g_EAAT2*eaat2_current(X(:,4),Glu_bath(t)', H_out, X(:,2), X(:,5),X(:,7)),'LineWidth',2); hold on;
ylabel('I_{EAAT2}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 tmaxplot])
    legend('Glu bath, calcium transient', ...
    'Glu bath, no calcium transient');

figure(12);
    if jj==1
        plot(t,g_Kir41.*kir41_current(X(:,5),X(:,2),X(:,7)),'LineWidth',4);
    elseif jj==2
        plot(t,g_Kir41.*kir41_current(X(:,5),X(:,2),X(:,7)),'--','LineWidth',4);
    else
        plot(t,g_Kir41.*kir41_current(X(:,5),X(:,2),X(:,7)),'-.','LineWidth',4);
    end
hold on
ylabel('I_{Kir}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 tmaxplot])
    legend('Glu bath, calcium transient', ...
    'Glu bath, no calcium transient');

figure(13);
if jj==1
        plot(t,g_Na_leak1.*na_leak_current(X(:,4),X(:,1), X(:,7)),'LineWidth',4);
    elseif jj==2
        plot(t,g_Na_leak1.*na_leak_current(X(:,4),X(:,1), X(:,7)),'--','LineWidth',4);
    else
        plot(t,g_Na_leak1.*na_leak_current(X(:,4),X(:,1), X(:,7)),'-.','LineWidth',4);
end
hold on
ylabel('I_{NaLeak}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 tmaxplot])
    legend('Glu bath, calcium transient', ...
    'Glu bath, no calcium transient');

% figure(13);
% V_NaCa = (R*T/F).*(3*log(X(:,4)./X(:,1))-log(X(:,6)./X(:,3))).*1e3; %mV
% plot(t,V_NaCa,'LineWidth',2); hold on;
% ylabel('V_{NaCa}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 tmaxplot])
% 
% figure(14);
% %assuming neuron has a steady concentration of [K+]i = 135 mM
% VKN = (R*T/F)*log(X(:,5)./135)*1e3; %mV
% plot(t,VKN,'LineWidth',2); hold on;
% ylabel('E_{K,N}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 tmaxplot])
% legend('Elevated potassium, calcium response', ...
%     'Elevated potassium, no calcium response',...
%     'Resting potassium');
% 
% figure(15);
% %assuming neuron has a steady concentration of [Na+]i = 12 mM
% VKN = (R*T/F)*log(X(:,4)./12)*1e3; %mV
% plot(t,VKN,'LineWidth',2); hold on;
% ylabel('E_{Na,N}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 tmaxplot])
% legend('Elevated potassium, calcium response', ...
%     'Elevated potassium, no calcium response',...
%     'Resting potassium');

Na_outs{jj} = X(:,4);
K_outs{jj} = X(:,5);
ts{jj} = t;

% save('Na_K_outs_Glubath.mat','ts','Na_outs', 'K_outs');
end
% figure(5);
