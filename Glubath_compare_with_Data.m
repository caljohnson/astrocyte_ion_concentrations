%Glubath_compare_with_Data.m
%compare glutamate bath stimulation with experimental data
%in particular the Na_in and Ca_in changes

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
g_EAAT2 = 5e-14; %tuned to get correct Nain and Cain responses as in Kirischuk et al 2007
g_Kir41 = 1e-14;
g_NCX = 1e-9;
rho_NKA = 1e-13;
VolE = 1.41e-18;  %L, from 1.41e-3 micrometers^3, extracellular/synaptic cleft volume from Handy, Lawley, Borisyuk 2018
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
K_out = 3.2; %rest in this system
Na_out = 140; %mM - Kirischuk et al. 2012
Ca_out = 0.9; %mM - rest state in this model
H_out = (10^(-7.4))*1e3; %mM - Kirischuk et al. 2012

Vm = -80; %mV, resting astrocyte membrane potential Kirischuk et al. 2012


Na_outs = {};
K_outs = {};
ts = {};
tmaxplot = 1e2;
tspan = [0,tmaxplot];

Glu_bath = @(t) 0*t + 1*(t<20); %1mM bath for 20 seconds

Ca_force = @(t) 0.*t;

   
X0 = [Na_in; K_in; Ca_in; Na_out; K_out; Ca_out; Vm];

options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',1e-3);
[t,X] = ode15s(@(t,X) system_eqns_with_Glu_bath(t,X,Glu_bath,Ca_force), tspan, X0,options);
dt = 1e-3;
t0 = tspan(1):dt:tspan(2);
X = interp1(t,X,t0);
t = t0;

labels = {'[Na^+]_i, mM'; '[K^+]_i, mM'; '[Ca^{2+}]_i, mM';...
    '[Na^+]_e, mM'; '[K^+]_e, mM'; '[Ca^{2+}]_e, mM';...
     'V_m, mV';};

%load data
[t_data1, Na_data] = csvimport('data/kirischuk_etal_2007_fig2a.csv',...
                    'columns', [1, 2] ,'noHeader', true); 
% t_data1 = t_data1-5; %shift by 5 seconds to match with stimulus start in exp
                
[t_data2, Ca_data] = csvimport('data/kirischuk_etal_2007_fig1a_100mM.csv',...
                    'columns', [1, 2] ,'noHeader', true); 
t_data2 = t_data2-5; %shift by 5 seconds to match with stimulus start in exp
Ca_data = Ca_data.*1e-6; %scale to nM
 
figure(1); clf
plot(t,X(:,1),'-','LineWidth',4); hold on
plot(t_data1, Na_data,'o','LineWidth',4);
ylabel(labels(1)); xlabel('t (sec)');
set(gca,'FontSize',20);
legend('model', ...
'Kirischuk et al. 2007 data');
xlim([0 tmaxplot])

figure(2); clf
plot(t,X(:,3),'-','LineWidth',4); hold on
plot(t_data2, Ca_data,'o','LineWidth',4);
ylabel(labels(3)); xlabel('t (sec)');
set(gca,'FontSize',20);
legend('model', ...
'Kirischuk et al. 2007 data');
xlim([0 tmaxplot])

figure(3); clf
plot(t,Glu_bath(t), 'LineWidth',4);
hold on
ylabel('Glu_{e}'); xlabel('t'); set(gca,'FontSize',20); xlim([0 tmaxplot])
                
