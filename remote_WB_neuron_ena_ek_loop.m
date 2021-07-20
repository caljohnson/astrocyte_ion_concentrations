%wang-buzsaki hippocampal neuron astrocyte effects timingsed current inject
%test to see if faster astrocytic potassium uptake and sodium dynamics have an effect of neural
%excitability

%demonstrates that elevated potassium leads to increased excitability
% and that Ca-transient helps lower excitability to compensate

addpath('./src'); %close all; clear;
% load('Na_K_outs.mat');
% load('Na_K_outs_superelevatedK.mat');
% load('Na_K_outs_onlyCatransient.mat');

F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature
% tmax = 1e3; %v1
tmax = 5e3; %v2,v3, v4
thresh = -5; %v4

%neural concentrations
K_in = 93.2; %mM  - set to make V_K = -90 mV at rest
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest

% %neural reversal potentials
% V_Na = (R*T/F).*log(Na_out./Na_in).*1e3;
% V_K = (R*T/F).*log(K_out./K_in).*1e3;

% E_Ks = -95:5:-50;
% E_Nas = 52.5:0.5:55.5;
% Iapps = 0:0.01:0.2; %v1
% Iapps = 0:0.005:0.11; %v2

%v3 params
E_Ks = -95:2.5:-50;
E_Nas = 52.5:0.25:55.5;
% Iapps = 0.10:0.01:0.12;
Iapps = 0.11; %v4

freq = zeros(size(Iapps,2),size(E_Ks,2),size(E_Nas,2));
delay = zeros(size(Iapps,2),size(E_Ks,2),size(E_Nas,2));
V0 = -64;
tic
for ii=1:size(Iapps,2)
    I = Iapps(ii);
    
    fprintf(['\r Iapp = ' num2str(I)]);
    
    for jj=1:size(E_Ks,2)
        E_K = E_Ks(jj);
        fprintf(['EK = ' num2str(E_K)]);
      for kk= 1:size(E_Nas,2)
            E_Na = E_Nas(kk);
            fprintf(['ENa = ' num2str(E_Na)]);

    X0 = [V0;0.78;0.088;];
    tvec = [0 tmax];
    [t,X] = ode23s(@(t,X) wang_buzsaki_hippocampal_neuron_ode(t,X,I,E_Na,E_K), tvec, X0);
    toc
%     figure(12); plot(t,X(:,1), 'LineWidth',4);
    spiked_1 = 0;
    spiked_2 = 0;
    tt = 1;
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
        delay(ii,jj,kk) = NaN;
    else
        delay(ii,jj,kk) = t(spike_time1);
    end
    
    if spiked_2 == 0
        freq(ii,jj,kk) = 0;
%         figure(ii);
%         plot(E_Na, E_K, 'rx'); hold on
    else
        freq(ii,jj,kk) = 1e3.*1/(t(spike_time2) - t(spike_time1));
%         figure(ii);
%         plot(E_Na, E_K, 'bo'); hold on
    end
      end
    end
%     figure(ii);
%     surf(E_Nas,E_Ks,freq); view(2); shading flat; colorbar
% xlabel('E_{Na} (mV)'); ylabel('E_K (mV)');
% title(['Frequency (Hz), I_{app} = ' num2str(Iapps(ii))]);
% set(gca,'FontSize',20);
% 
%     figure(ii+3);
%     surf(E_Nas,E_Ks,delay); view(2); shading flat; colorbar
% xlabel('E_{Na} (mV)'); ylabel('E_K (mV)');
% title(['Delay (msec), I_{app} = ' num2str(Iapps(ii))]);
% set(gca,'FontSize',20);

 


end

% save('remote_WB_neuron_ena_ek_loop.mat'); %v1
% save('remote_WB_neuron_ena_ek_loop_v2.mat'); 
% save('remote_WB_neuron_ena_ek_loop_v3.mat');
save('remote_WB_neuron_ena_ek_loop_v4.mat');