%wang-buzsaki hippocampal neuron astrocyte effects timingsed current inject
%test to see if faster astrocytic potassium uptake and sodium dynamics have an effect of neural
%excitability

%demonstrates that elevated potassium leads to increased excitability
% and that Ca-transient helps lower excitability to compensate

addpath('./src'); close all; clear;
% load('Na_K_outs.mat');
% load('Na_K_outs_superelevatedK.mat');
% load('Na_K_outs_onlyCatransient.mat');

F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature
tmax = 1e2;

%neural concentrations
K_in = 93.2; %mM  - set to make V_K = -90 mV at rest
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest

% %neural reversal potentials
% V_Na = (R*T/F).*log(Na_out./Na_in).*1e3;
% V_K = (R*T/F).*log(K_out./K_in).*1e3;

E_Ks = -95:1:-50;
E_Nas = 52.5:0.5:55.5;
% Iapps = 0.12:0.005:0.13;
Iapps = 0.11;

thresh = -5; %firing threshold

[x,y] = meshgrid(E_Nas,E_Ks);
x = x(:);
y = y(:);

freq = zeros(size(E_Ks,2),size(E_Nas,2));
delay = zeros(size(E_Ks,2),size(E_Nas,2));
% V0 = -64;
V0 = -50;
for ii=1:size(Iapps,2)
    I = Iapps(ii);
    for jj=1:size(E_Ks,2)
        E_K = E_Ks(jj);
      for kk= 1:size(E_Nas,2)
            E_Na = E_Nas(kk);

    X0 = [V0;0.78;0.088;];
    tvec = [0 tmax];
    tic
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
        delay(jj,kk) = NaN;
    else
        delay(jj,kk) = t(spike_time1);
    end
    
    if spiked_2 == 0
        freq(jj,kk) = 0;
%         figure(ii);
%         plot(E_Na, E_K, 'rx'); hold on
    else
        freq(jj,kk) = 1e3.*1/(t(spike_time2) - t(spike_time1));
%         figure(ii);
%         plot(E_Na, E_K, 'bo'); hold on
    end
      end
    end
    figure(ii);
    surf(E_Nas,E_Ks,freq); view(2); shading flat; 
%     scatter(x,y,30,freq(:),'LineWidth',4);
    colorbar
    xlabel('E_{Na} (mV)'); ylabel('E_K (mV)');
    xlim([min(E_Nas) max(E_Nas)]); ylim([min(E_Ks) max(E_Ks)]);
    title(['Frequency (Hz), I_{app} = ' num2str(Iapps(ii))]);
    set(gca,'FontSize',20);

    figure(ii+3);
    surf(E_Nas,E_Ks,delay); view(2); shading flat;
%     scatter(x,y,30,delay(:),'LineWidth',4);
    colorbar
    xlabel('E_{Na} (mV)'); ylabel('E_K (mV)');
    xlim([min(E_Nas) max(E_Nas)]); ylim([min(E_Ks) max(E_Ks)]);
    title(['Delay (msec), I_{app} = ' num2str(Iapps(ii))]);
    set(gca,'FontSize',20);

end

