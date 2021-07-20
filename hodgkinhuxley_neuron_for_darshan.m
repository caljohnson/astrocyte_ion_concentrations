%Hodgkin-Huxley type neuron code demo

addpath('./src'); close all; clear;


%parameters
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature
V_Na = 55; %mV
V_K = -90; %mV
V0 = -60; %mV - initial condition
X0 = [V0;1;0;]; %full set of initial conditions
 
%ode timespan
t0 = 0;
tmax = 1e2;
tvec = [t0 tmax];

%Applied current - can make this whatever you want to have an "input"
%signal
I = @(t) 0.*t;% + 0.2.*(t>50 && t<100)+10.*(t>250 && t<300);

%run ode23s in this case because the neuron oscillates
[t,X] = ode23s(@(t,X) wang_buzsaki_hippocampal_neuron_ode(t,X,I,V_Na,V_K), tvec, X0);    
figure(1); hold on
plot(t,X(:,1),'LineWidth',4); ylim([-100,40])   
set(gca,'FontSize',20);
xlabel('time (seconds)'); ylabel('neuron membrane potential (mV)');


V0 = -58; %mV - initial condition that gets it oscillating
X0 = [V0;1;0;]; %full set of initial conditions
[t,X] = ode23s(@(t,X) wang_buzsaki_hippocampal_neuron_ode(t,X,I,V_Na,V_K), tvec, X0);    
figure(1); plot(t,X(:,1),'LineWidth',4); 



