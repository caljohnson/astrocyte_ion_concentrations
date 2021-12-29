%wang-buzsaki hippocampal neuron base-case test
%test to check whether XPPAUT and matlab give different answers

addpath('./src');
% close all; clear;

%params
V_K = -90; %eq potential
V_Na = 55; %eq potential
% V_K = -58;
% V_Na = 55;
% I_app = 0.165; %applied current
I_app = 0.5;

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

spiked_1 = 0;
spiked_2 = 0;
thresh= 0;
tt = round(size(t,1)/2); %start looking for frequency in second half of simulation trace
while tt <= size(t,1) && spiked_2 == 0
    if spiked_1 == 1 && X(tt,1) > thresh && X(tt-1,1) <=thresh && spiked_2 == 0
        spiked_2 = 1;
        spike_time2 = tt;
    end
    if X(tt,1) > thresh && X(tt-1,1) <=thresh && spiked_1 == 0
        spiked_1 = 1;
        spike_time1 = tt;
    end
    tt = tt+1;
end

if spiked_2 == 0
    freq = 0
else
    freq = 1e3.*1/(t(spike_time2) - t(spike_time1))
end


