%wang-buzsaki neuron E_Na, E_K trajectory simulation

addpath('./src'); close all; clear;
load('Na_K_outs_Glupulses.mat');

F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature
tmax = 5e4;

%neural concentrations
K_in = 93.2; %mM  - set to make V_K = -90 mV at rest
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest

%external concentrations
K_out = K_outs{1};
Na_out = Na_outs{1};

% I = 0.165;

I = 0.16

tsmax = ts{1}(end);
Na_e = @(t) interp1(ts{1},Na_out, t*1e-3,'linear','extrap').*(t*1e-3<=tsmax) + Na_out(end).*(t*1e-3>tsmax);
K_e = @(t) interp1(ts{1},K_out, t*1e-3,'linear','extrap').*(t*1e-3<=tsmax) + K_out(end).*(t*1e-3>tsmax);

E_Na = @(t) (R*T/F).*log(Na_e(t)./Na_in).*1e3;
E_K = @(t) (R*T/F).*log(K_e(t)./K_in).*1e3;

thresh = 0; %firing threshold

V0 = -64;

X0 = [V0;0.78;0.088;];
tvec = [0 tmax];
tic
% options = odeset('MaxStep',0.1);
options = odeset();
[t,X] = ode23s(@(t,X) wb_neuron_ode_variable_es(t,X,I,E_Na,E_K), tvec, X0,options);
toc

figure(1); clf
subplot(3,1,1);
plot(t.*1e-3,X(:,1),'LineWidth',4); hold on
xlabel('time, sec'); ylabel('V_n, mV');
set(gca, 'FontSize',20);

subplot(3,1,2);
plot(t.*1e-3,Na_e(t),'LineWidth',4); hold on
xlabel('time, sec'); ylabel('[Na^+]_e, mM');
set(gca, 'FontSize',20);

subplot(3,1,3);
plot(t.*1e-3,K_e(t),'LineWidth',4); hold on
xlabel('time, sec'); ylabel('[K^+]_e, mM');
set(gca, 'FontSize',20);

return
tic
[t2,X2] = ode23s(@(t,X) wang_buzsaki_hippocampal_neuron_ode(t,X,I,E_Na(tsmax),E_K(tsmax)), tvec, X0,options);
toc

subplot(3,1,1);
plot(t2.*1e-3,X2(:,1),':','LineWidth',4);

subplot(3,1,2);
plot(t2.*1e-3,Na_e(t(end))+0.*t2,':','LineWidth',4);

subplot(3,1,3);
plot(t2.*1e-3,K_e(t(end))+0.*t2,':','LineWidth',4);