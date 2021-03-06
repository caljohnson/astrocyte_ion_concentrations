%wang-buzsaki hippocampal neuron astrocyte effects timingsed current inject
%test to see if faster astrocytic potassium uptake and sodium dynamics have an effect of neural
%excitability

%demonstrates that elevated potassium leads to increased excitability
% and that Ca-transient helps lower excitability to compensate

% load('remote_WB_neuron_ena_ek_loop_v1.mat'); %v1 - EK-95 to -50, ENa50 to 60
% load('remote_WB_neuron_ena_ek_loop_v2.mat'); %v2 - -95 to -80, ENa50 to 60
load('remote_WB_neuron_ena_ek_loop_v3.mat'); %v2 - EK-95 to -80, ENa50 to 60, Iapp 0.1-0.2

[x,y] = meshgrid(E_Nas,E_Ks);
x = x(:);
y = y(:);

for ii=1:size(Iapps,2)
         
    freqs = reshape(freq(ii,:,:),size(E_Ks,2),size(E_Nas,2));
    delays = reshape(delay(ii,:,:),size(E_Ks,2),size(E_Nas,2));
   
    freqs(freqs==0) = NaN;
    figure(ii);
%     surf(E_Nas,E_Ks,freqs); view(2); shading flat; 
    scatter(x,y,30,freqs(:),'LineWidth',4); hold on
    colorbar
xlabel('E_{Na} (mV)'); ylabel('E_K (mV)');
xlim([min(E_Nas) max(E_Nas)]); ylim([min(E_Ks) max(E_Ks)]);
title(['Frequency (Hz), I_{app} = ' num2str(Iapps(ii))]);
set(gca,'FontSize',20);

%     figure(ii+size(Iapps,2));
% %     surf(E_Nas,E_Ks,delays); view(2); shading flat; colorbar
%     scatter(x,y,30,delays(:),'LineWidth',4);
%     colorbar
% xlabel('E_{Na} (mV)'); ylabel('E_K (mV)');
% xlim([min(E_Nas) max(E_Nas)]); ylim([min(E_Ks) max(E_Ks)]);
% title(['Delay (msec), I_{app} = ' num2str(Iapps(ii))]);
% set(gca,'FontSize',20);

 
end
% 
% % for ii=1:2*size(Iapps,2)
% for ii=1:size(Iapps,2)
%     figure(ii); set(gcf, 'Position', get(0, 'Screensize'));
% end