%figs_WB_neuron_ena_ek_astrocyte_trajectory

load('remote_WB_neuron_ena_ek_loop_v4.mat');
% load('remote_WB_neuron_ena_ek_loop_v5.mat');
% load('remote_WB_neuron_ena_ek_loop_v6.mat');
% load('remote_WB_neuron_ena_ek_loop_v3.mat');

close all; 
[x,y] = meshgrid(E_Nas,E_Ks);
x = x(:);
y = y(:);

ii = 8; %corresponds to Iapp =0.17

freqs = reshape(freq(ii,:,:),size(E_Ks,2),size(E_Nas,2));
freqs(freqs==0)= NaN;
delays = reshape(delay(ii,:,:),size(E_Ks,2),size(E_Nas,2));
   
figure(1);
% surf(E_Nas,E_Ks,freqs); view(2); shading flat; 
scatter(x,y,100,freqs(:),'LineWidth',4);
hold on; colorbar
% caxis([30 100]);
xlabel('E_{Na} (mV)'); ylabel('E_K (mV)');
xlim([min(E_Nas) max(E_Nas)]); ylim([min(E_Ks) max(E_Ks)]);
title(['Frequency (Hz), I_{app} = ' num2str(Iapps(ii))]);
set(gca,'FontSize',20); set(gcf, 'Position', get(0, 'Screensize'));

% figure(2);
% % surf(E_Nas,E_Ks,delays); view(2); shading flat;
% scatter(x,y,100,delays(:),'LineWidth',4); 
% hold on; colorbar
% xlabel('E_{Na} (mV)'); ylabel('E_K (mV)');
% xlim([min(E_Nas) max(E_Nas)]); ylim([min(E_Ks) max(E_Ks)]);
% title(['Delay (msec), I_{app} = ' num2str(Iapps(ii))]);
% set(gca,'FontSize',20); set(gcf, 'Position', get(0, 'Screensize'));

%highlight "normal" conditions
figure(1);
plot(55.0662, -90.0758, 'go', 'LineWidth',4,'MarkerSize',20);
%highlight elevated K conditions
% plot(55.0572, -54.7297, 'ro', 'LineWidth',4,'MarkerSize',20);

%highlight trajectories "Glu puff to rest"
% load('Na_K_outs_Glupulses.mat');
load('Na_K_outs_Glupulses.mat');
% load('Na_K_outs_2Glupulses.mat');
% load('Na_K_outs_Glu_diffpulses.mat');
K_out = K_outs{1};
Na_out = Na_outs{1};
K_in = 93.2; %mM  - set to make V_K = -90 mV at rest
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest
%  t = ts{1};
V_Nas = (R*T/F).*log(Na_out./Na_in).*1e3;
V_Ks = (R*T/F).*log(K_out./K_in).*1e3;
plot(V_Nas,V_Ks,'b-','LineWidth',4)

K_out = K_outs{2};
Na_out = Na_outs{2};
%  t = ts{1};
V_Nas = (R*T/F).*log(Na_out./Na_in).*1e3;
V_Ks = (R*T/F).*log(K_out./K_in).*1e3;
plot(V_Nas,V_Ks,'r-','LineWidth',4)

%with glu bath
% h= legend('nonzero frequencies due to step current (no Glu)', ...
%     'step current (no Glu), resting conditions',...
%     'step current as Glu bath, astrocyte effects w/ Calcium transient',...
%     'step current as Glu bath, astrocyte effects w/o Calcium transient');

% h= legend('nonzero frequencies due to step current (no Glu)', ...
%     'step current (no Glu), resting conditions',...
%     'initial 3mM Glu, astrocyte effects w/ Calcium transient',...
%     'initial 3mM Glu, astrocyte effects w/o Calcium transient');
% set(h,'Location','northwest');


% load('Na_K_outs_Glupulses.mat');
load('Na_K_outs_2Glupulses.mat');
% load('Na_K_outs_Glu_diffpulses.mat');
K_out = K_outs{1};
Na_out = Na_outs{1};
K_in = 93.2; %mM  - set to make V_K = -90 mV at rest
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest
%  t = ts{1};
V_Nas = (R*T/F).*log(Na_out./Na_in).*1e3;
V_Ks = (R*T/F).*log(K_out./K_in).*1e3;
plot(V_Nas,V_Ks,'b--','LineWidth',4)

K_out = K_outs{2};
Na_out = Na_outs{2};
%  t = ts{1};
V_Nas = (R*T/F).*log(Na_out./Na_in).*1e3;
V_Ks = (R*T/F).*log(K_out./K_in).*1e3;
plot(V_Nas,V_Ks,'r--','LineWidth',4)

%with glu bath
% h= legend('nonzero frequencies due to step current (no Glu)', ...
%     'step current (no Glu), resting conditions',...
%     'step current as Glu bath, astrocyte effects w/ Calcium transient',...
%     'step current as Glu bath, astrocyte effects w/o Calcium transient');

% h= legend('nonzero frequencies due to step current (no Glu)', ...
%     'step current (no Glu), resting conditions',...
%     'initial 3mM Glu, astrocyte effects w/ Calcium transient',...
%     'initial 3mM Glu, astrocyte effects w/o Calcium transient',...
%     '2 pulses 3mM Glu, astrocyte effects w/ Calcium transient',...
%     '2 pulses 3mM Glu, astrocyte effects w/o Calcium transient');
% set(h,'Location','northwest');

%3 pulses

load('Na_K_outs_3Glupulses.mat');
K_out = K_outs{1};
Na_out = Na_outs{1};
K_in = 93.2; %mM  - set to make V_K = -90 mV at rest
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest
%  t = ts{1};
V_Nas = (R*T/F).*log(Na_out./Na_in).*1e3;
V_Ks = (R*T/F).*log(K_out./K_in).*1e3;
plot(V_Nas,V_Ks,'b:','LineWidth',4)

K_out = K_outs{2};
Na_out = Na_outs{2};
%  t = ts{1};
V_Nas = (R*T/F).*log(Na_out./Na_in).*1e3;
V_Ks = (R*T/F).*log(K_out./K_in).*1e3;
plot(V_Nas,V_Ks,'r:','LineWidth',4)

h= legend('nonzero frequencies due to step current (no Glu)', ...
    'step current (no Glu), resting conditions',...
    'initial 3mM Glu, astrocyte effects w/ Calcium transient',...
    'initial 3mM Glu, astrocyte effects w/o Calcium transient',...
    '2 pulses 3mM Glu, astrocyte effects w/ Calcium transient',...
    '2 pulses 3mM Glu, astrocyte effects w/o Calcium transient',...
    '3 pulses 3mM Glu, astrocyte effects w/ Calcium transient',...
    '3 pulses 3mM Glu, astrocyte effects w/o Calcium transient');
set(h,'Location','northwest');

return
%highlight trajectories "elevated K to rest"
load('Na_K_outs_superelevatedK.mat')
%condition 1 = elevated potassium, calcium transient
K_out = K_outs{1};
Na_out = Na_outs{1};
K_in = 93.2; %mM  - set to make V_K = -90 mV at rest
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest
%  t = ts{1};
V_Nas = (R*T/F).*log(Na_out./Na_in).*1e3;
V_Ks = (R*T/F).*log(K_out./K_in).*1e3;
plot(V_Nas,V_Ks,'m-','LineWidth',4)
% plot_dir(V_Nas,V_Ks,'m')

%condition 2 = elevated potassium, no Ca-transient
K_out = K_outs{2};
Na_out = Na_outs{2};
%  t = ts{2};
V_Nas = (R*T/F).*log(Na_out./Na_in).*1e3;
V_Ks = (R*T/F).*log(K_out./K_in).*1e3;
plot(V_Nas,V_Ks,'c-','LineWidth',4)
plot_dir(V_Nas,V_Ks,'c')

h= legend('nonzero frequencies', 'resting conditions','elevated K^+',...
    'astrocyte effects, w/ Calcium transient',...
    'astrocyte effects, no Calcium transient');
set(h,'Location','northwest');


