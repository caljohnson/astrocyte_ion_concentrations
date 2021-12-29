%compare WB and EIF neuron models
% with the "same" astrocyte effects
% at different constant applied current levels

%EIF neuron test

clear; clc; close all

addpath('./');

global V_T del_T g_L g_T E_L C

I_app = 0.165;

%wang-buzsaki hippocampal neuron - astrocyte effects
load('Na_K_outs_Glupulses.mat');

%neural concentrations
K_in = 93.2; %mM  - set to make V_K = -90 mV at rest
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest

%external concentrations
K_out = K_outs{1};
Na_out = Na_outs{1};

%more params
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

tsmax = ts{1}(end);
Na_e = @(t) interp1(ts{1},Na_out, t*1e-3,'linear','extrap').*(t*1e-3<=tsmax) + Na_out(end).*(t*1e-3>tsmax);
K_e = @(t) interp1(ts{1},K_out, t*1e-3,'linear','extrap').*(t*1e-3<=tsmax) + K_out(end).*(t*1e-3>tsmax);

E_Na = @(t) (R*T/F).*log(Na_e(t)./Na_in).*1e3;
E_K = @(t) (R*T/F).*log(K_e(t)./K_in).*1e3;

%initial condition
V0 = -55;
X0 = [V0;0.78;0.088;];
    
%solve and plot    
tvec = [0 5e4];
tic
[t,X] = ode23s(@(t,X) wb_neuron_ode_variable_es(t,X,I_app,E_Na,E_K), tvec, X0);
toc

figure(1); hold on
plot(t.*1e-3,X(:,1),'-','LineWidth',4); ylim([-100,40])   
set(gca,'FontSize',20);
xlabel('time (sec)'); ylabel('neuron membrane potential (mV)');

%wang-buzsaki hippocampal neuron - current effects
%params
V_K = -90; %eq potential
V_Na = 55; %eq potential

%total current - applied + astrocyte-mediated current
%muA/cm^2
% I_ast = @(t) E_Na(t) - 55.0662;
% I_ast = @(t) 0.*t -0.5e-3.*(t-4e3).*(t>4e3).*(t<6e3)...
%     -1.*(t>6e3).*(t<14e3)+ (-1-0.333e-3.*(14e3-t)).*(t>14e3).*(t<17e3);
I_ast = @(t) 0.*t -0.5e-3.*(t-4e3).*(t>4e3).*(t<6e3)...
    -1.*(t>6e3).*(t<13.5e3)+ (-1-0.25e-3.*(13.5e3-t)).*(t>13.5e3).*(t<15.5e3)...
    -0.5.*(t>15.5e3).*(t<18e3)+(-0.5-0.25e-3.*(18e3-t)).*(t>18e3).*(t<20e3);
I = @(t) 0.*t + I_app + 0.005.*I_ast(t);

%initial condition
V0 = -55;
X0 = [V0;0.78;0.088;];
    
%solve and plot    
tvec = [0 5e4];
tic
[t,X] = ode23s(@(t,X) WB_neuron_ode_variable_Icurrent(t,X,I,V_Na,V_K), tvec, X0);
toc

figure(1); hold on
plot(t.*1e-3,X(:,1),'--','LineWidth',4); ylim([-100,40])   
set(gca,'FontSize',20);
xlabel('time (sec)'); ylabel('neuron membrane potential (mV)');



figure(2);
subplot(3,1,1)
plot(t.*1e-3,E_Na(t),'-','LineWidth',4);
xlabel('time (sec)'); ylabel('E_{Na} (mV)');
set(gca,'FontSize',20);
subplot(3,1,2)
plot(t.*1e-3,E_K(t),'-','LineWidth',4);
xlabel('time (sec)'); ylabel('E_{K} (mV)');
set(gca,'FontSize',20);
subplot(3,1,3)
plot(t.*1e-3,I(t),'-','LineWidth',4); hold on
plot(t.*1e-3, I_app + 0.003.*I_ast(t),'--','LineWidth',4);
set(gca,'FontSize',20);
xlabel('time (sec)'); ylabel('I_{app}+I_{ast} (\muA/cm^2)');
legend('WB','EIF');


%EIF neuron - paramaters from Fourcaud et al. 2003
C = 1; %muF/cm^2, membrane capacitance
g_L = 0.1;%mS/cm^2, leakage current conductance
g_T = g_L;%threshold current conductance 
del_T = 3.48; %mV, slope factor 
E_L = -65; %mV, resting potential
V_T = -59.9; %mV, soft threshold potential
V_r = -68; %mV, reset potential

%astrocyte-mediated current, reduces excitability
%muA/cm^2
% I_ast = @(t) (E_Na(t) - 55.0662).*3e-3;
% I_ast = @(t) 0.*t -0.5e-3.*(t-4e3).*(t>4e3).*(t<6e3)...
%     -1.*(t>6e3).*(t<14e3)+ (-1-0.5e-3.*(14e3-t)).*(t>14e3).*(t<16e3);

%total current - applied + astrocyte-mediated current
%muA/cm^2
I = @(t) 0.*t + I_app + 0.003.*I_ast(t);

t0 = 0;
tmax = 1e3;
tvec = [t0 tmax];
tmax_final = 5e4;


X0 = V_r;

tic
options = odeset('Events',@depolEvents);
warning('off','all')
clear t X
for k= 1:200
    [t{k},X{k}] = ode23s(@(t,X) EIF_neuron_ode(t,X,I), tvec, X0,options);
    if t{k}(end) == tvec(end)
        X0 = X{k}(end);
    else
        X0 = V_r;
    end
    if t{k}(end) > tmax_final
        break
    else
        tvec = [t{k}(end) t{k}(end)+tmax];
    end
end
toc
tv = cat(1,t{:});
xv = cat(1,X{:});

figure(1); hold on
plot(tv.*1e-3,xv(:,1),'-.','LineWidth',4); ylim([-100,40])
set(gca,'FontSize',20);
xlabel('time (sec)'); ylabel('neuron membrane potential (mV)');


legend('WB with (astrocyte-caused) E_{Na}, E_{K} changes',...
    'WB with astrocyte-current',...
    'EIF with astrocyte-current');
xlim([0 30]);


function [position,isterminal,direction] = depolEvents(t,y)
global V_T
    if isnan(y)
        position = 0;
    elseif isinf(y)
        position = 0;
    elseif y>=V_T
        position = 0;
    else
        position = 1; % The value that we want to be zero cutoff
    end
    isterminal = 1;  % Halt integration 
    direction = 1;   % The zero can be approached from either direction
end
