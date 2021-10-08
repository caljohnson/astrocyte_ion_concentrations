%wang-buzsaki hippocampal neuron freq
%how does applied current influence frequency?

addpath('./src'); close all; clear;


F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

V_K = -90; %resting
V_Na = 55; %resting

%simulation times
tmax = 1e2;
tvec = [0 tmax];
tic
Iapps = linspace(-0.7,-0.5,200); %applied current values to loop over
X0 = [-55;0.78;0.088;]; %initial condition
for kk=1:size(Iapps,2)
    
    I = Iapps(kk);
    
    [t,X] = ode23s(@(t,X) wang_buzsaki_hippocampal_neuron_ode(t,X,I,V_Na,V_K), tvec, X0);
%     figure(kk); hold on
%     plot(t,X(:,1),'LineWidth',4); ylim([-100,40])   
%     set(gca,'FontSize',20);
%     xlabel('time (msec)'); ylabel('neuron membrane potential (mV)');
   
    spiked_1 = 0;
    spiked_2 = 0;
    tt = 1;
    thresh = -20;
    while tt <= size(t,1) && spiked_2 == 0
        if spiked_1 == 1 && X(tt,1) > thresh && X(tt-1) <=thresh && spiked_2 == 0
            spiked_2 = 1;
            spike_time2 = tt;
        end
        if X(tt,1) > thresh && X(tt-1) <=thresh && spiked_1 == 0
            spiked_1 = 1;
            spike_time1 = tt;
        end
        tt = tt+1;
    end

    if spiked_1 == 0
        delay = NaN;
    else
        delay = t(spike_time1);
    end

    if spiked_2 == 0
        freq = 0;
    else
        freq = 1e3.*1/(t(spike_time2) - t(spike_time1));
    end

    freqs(kk) = freq;
end
toc
% figure(kk+1);
figure(1);
plot(Iapps, freqs, 'o', 'LineWidth',4)
xlabel('I_{app} (mA/cm^2)'); ylabel('Frequency (Hz)');

