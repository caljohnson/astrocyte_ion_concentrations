%figs_WB_neuron_ena_ek_astrocyte_trajectory

load('remote_WB_neuron_ena_ek_loop_v4.mat');

close all; 
[x,y] = meshgrid(E_Nas,E_Ks);
x = x(:);
y = y(:);
         
for ii=1:size(Iapps,2)
freqs = reshape(freq(ii,:,:),size(E_Ks,2),size(E_Nas,2));
freqs(freqs==0)= NaN;
delays = reshape(delay(ii,:,:),size(E_Ks,2),size(E_Nas,2));
   
figure(ii);
% surf(E_Nas,E_Ks,freqs); view(2); shading flat; 
scatter(x,y,100,freqs(:),'LineWidth',4);
hold on; colorbar
xlabel('E_{Na} (mV)'); ylabel('E_K (mV)');
xlim([min(E_Nas) max(E_Nas)]); ylim([min(E_Ks) max(E_Ks)]);
title(['Frequency (Hz), I_{app} = ' num2str(Iapps(ii))]);
set(gca,'FontSize',20); set(gcf, 'Position', get(0, 'Screensize'));

%highlight "normal" conditions
plot(55.0662, -90.0758, 'go', 'LineWidth',4,'MarkerSize',20);
%highlight elevated K conditions
plot(55.0572, -54.7297, 'ro', 'LineWidth',4,'MarkerSize',20);

%highlight trajectories "elevated K to rest"
load('Na_K_outs_elevatedK_noCa.mat')

%condition 1 = elevated potassium
K_out = K_outs{1};
Na_out = Na_outs{1};
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest
%  t = ts{1};
V_Nas = (R*T/F).*log(Na_out./Na_in).*1e3;
V_Ks = (R*T/F).*log(K_out./K_in).*1e3;
plot(V_Nas,V_Ks,'m-','LineWidth',4)
% plot_dir(V_Nas,V_Ks,'m')


h= legend('nonzero frequencies', 'resting conditions','elevated K^+',...
...%     'astrocyte effects, w/ Calcium transient',...
    'astrocyte modulatory trajectory');
set(h,'Location','northwest');

end

