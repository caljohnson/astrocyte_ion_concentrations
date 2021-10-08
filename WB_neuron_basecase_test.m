%wang-buzsaki hippocampal neuron base-case test
%test to check whether XPPAUT and matlab give different answers

addpath('./src');
% close all; clear;

%params
% V_K = -90; %eq potential
% V_Na = 55; %eq potential
V_K = -58;
V_Na = 55;
I_app = 0.165; %applied current

%initial condition
V0 = -55;
X0 = [V0;0.78;0.088;];
    
%solve and plot    
tvec = [0 1e3];
[t,X] = ode23s(@(t,X) wang_buzsaki_hippocampal_neuron_ode(t,X,I_app,V_Na,V_K), tvec, X0);

figure(1); hold on
plot(t,X(:,1),'-','LineWidth',4); ylim([-100,40])   
set(gca,'FontSize',20);
xlabel('time (msec)'); ylabel('neuron membrane potential (mV)');


