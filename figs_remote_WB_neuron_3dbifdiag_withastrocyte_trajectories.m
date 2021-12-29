%wang-buzsaki hippocampal neuron astrocyte effects timingsed current inject
%test to see if faster astrocytic potassium uptake and sodium dynamics have an effect of neural
%excitability

%demonstrates that elevated potassium leads to increased excitability
% and that Ca-transient helps lower excitability to compensate


clear; close all; clc;
% load('remote_WB_neuron_ena_ek_loop_v4.mat'); %v4 - EK-95 to -50, ENa50 to 60, Iapp 0.1-0.2
%v5 has fixed frequency-finding algorithm
load('remote_WB_neuron_ena_ek_loop_v5.mat'); %v5 - EK-95 to -50, ENa50 to 60, Iapp 0.1-0.2 

load('Na_K_outs_Glupulses.mat');
K_out = K_outs{1};
Na_out = Na_outs{1};
K_in = 93.2; %mM  - set to make V_K = -90 mV at rest
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest
V_Nas = (R*T/F).*log(Na_out./Na_in).*1e3;
V_Ks = (R*T/F).*log(K_out./K_in).*1e3;

K_out2 = K_outs{2};
Na_out2 = Na_outs{2};
V_Nas2 = (R*T/F).*log(Na_out2./Na_in).*1e3;
V_Ks2 = (R*T/F).*log(K_out2./K_in).*1e3;

[x,y] = meshgrid(E_Nas,E_Ks);
x = x(:);
y = y(:);

for ii=1:size(Iapps,2)
         
    freqs = reshape(freq(ii,:,:),size(E_Ks,2),size(E_Nas,2));
    delays = reshape(delay(ii,:,:),size(E_Ks,2),size(E_Nas,2));
   
    freqs(freqs==0) = NaN;
    figure(ii); clf;
%     surf(E_Nas,E_Ks,freqs); view(2); shading flat; 
    scatter(x,y,30,freqs(:),'LineWidth',4); hold on
    colorbar
xlabel('E_{Na} (mV)'); ylabel('E_K (mV)');
xlim([min(E_Nas) max(E_Nas)]); ylim([min(E_Ks) max(E_Ks)]);
title(['Frequency (Hz), I_{app} = ' num2str(Iapps(ii))]);
set(gca,'FontSize',20);

%find abrupt changes in frequency (bifurcation line)
freqs(isnan(freqs)) = -1e10;
%logical indexes of the bifurcation line
bif_inds = ischange(freqs,2,'Threshold',1e10);
bif_inds2 = ischange(freqs,1,'Threshold',1e10);
bif_inds = bif_inds | bif_inds2;

% bif_inds = ischange(freqs,2,'MaxNumChanges',2);
x_bif{ii} = x(bif_inds);
y_bif{ii} = y(bif_inds);
I_bif{ii} = Iapps(ii).*ones(size(x(bif_inds)));
%plot x-y coords of the bifurcation line (roughly)
figure(ii);
% plot(x(bif_inds)-0.1250,y(bif_inds)-0.1250,'-','LineWidth',4);
% plot(x_bif{ii},y_bif{ii},'ro','LineWidth',4);

%get other bifurcation line
% bif_inds2 = ischange(freqs,1,'Threshold',500);
% bif_inds2 = ischange(freqs,1,'MaxNumChanges',1);
% x_bif2{ii} = x(bif_inds2);
% y_bif2{ii} = y(bif_inds2);
% I_bif2{ii} = Iapps(ii).*ones(size(x(bif_inds2)));
% plot(x_bif2{ii},y_bif2{ii},'mo','LineWidth',4);
%  
% %add to 3D bifurcation diagram
% figure(30);
% % scatter3(x(bif_inds)-0.1250, ...
% %     y(bif_inds)-0.1250,...
% %     Iapps(ii).*ones(size(x(bif_inds))),'b'); hold on
% scatter3(x(bif_inds), ...
%     y(bif_inds),...
%     Iapps(ii).*ones(size(x(bif_inds))),'b'); hold on
% xlabel('E_{Na} (mV)'); ylabel('E_K (mV)'); zlabel('I_{app} (mA)');
% zlim([0,1]);
% 


plot(V_Nas,V_Ks,'b-','LineWidth',4)
plot(V_Nas2,V_Ks2,'r-','LineWidth',4)



end

%add to 3D bifurcation diagram
figure(30);
x = [];
y = [];
z = [];
for ii=1:size(Iapps,2)
    scatter3(x_bif{ii},...
    y_bif{ii},...
    I_bif{ii},'b',...
    'HandleVisibility','off'); hold on

    xlabel('E_{Na} (mV)'); ylabel('E_K (mV)'); zlabel('I_{app} (mA)');
    zlim([0.1,0.3]);
    
    x = [x; x_bif{ii}];
    y = [y; y_bif{ii}];
    z = [z; I_bif{ii}];

end

% EKmax = -56;
EKmax = -60.5;
x1 = x(y<EKmax);
y1 = y(y<EKmax);
z1 = z(y<EKmax);

EKmax = -60.5;
zmax = 0.122;
x2 = x(y>=EKmax & z>zmax);
y2 = y(y>=EKmax & z>zmax);
z2 = z(y>=EKmax & z>zmax);


nx = 100 ; ny = 100 ;
xi = linspace(min(x1),max(x1),nx) ; 
yi = linspace(min(y1),max(y1),ny) ; 
[X,Y] = meshgrid(xi,yi) ;
F = scatteredInterpolant(x1,y1,z1) ; 
Z = F(X,Y);
surf(X,Y,Z, 'FaceColor','g', 'FaceAlpha',0.5, 'EdgeColor','none'); hold on


nx = 100 ; ny = 100 ;
xi = linspace(min(x2),max(x2),nx) ; 
yi = linspace(min(y2),max(y2),ny) ; 
[X2,Y2] = meshgrid(xi,yi) ;
F2 = scatteredInterpolant(x2,y2,z2) ;
Z2 = F2(X2,Y2) ;
surf(X2,Y2,Z2, 'FaceColor','r', 'FaceAlpha',0.5, 'EdgeColor','none')
zlim([0.1, 0.2]);
set(gca,'FontSize',20);

%add astrocyte trajectory at bottom of 3d plot
plot3(V_Nas,V_Ks,0.1.*ones(size(V_Nas)),'b:','LineWidth',4);
plot3(V_Nas2,V_Ks2,0.1.*ones(size(V_Nas2)),'r-','LineWidth',4);

legend('SNIC','HB',...
'Astrocyte modulatory trajectory, with Ca-transient',...
'Astrocyte modulatory trajectory, no Ca-transient',...
'Location','northeast');

%highlight trajectories "elevated K to rest"
load('Na_K_outs_elevatedK.mat')
%condition 1 = elevated potassium, calcium transient
K_out = K_outs{1};
Na_out = Na_outs{1};
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest
K_in = 93.2; %mM  - set to make V_K = -90 mV at restc
F = 96485; %C/mol, Faraday's constant
V_Nas = (R*T/F).*log(Na_out./Na_in).*1e3;
V_Ks = (R*T/F).*log(K_out./K_in).*1e3;
plot3(V_Nas,V_Ks,0.1.*ones(size(V_Nas)),'m-','LineWidth',4)

%condition 2 = elevated potassium, no Ca-transient
K_out = K_outs{2};
Na_out = Na_outs{2};
V_Nas2 = (R*T/F).*log(Na_out./Na_in).*1e3;
V_Ks2 = (R*T/F).*log(K_out./K_in).*1e3;
plot3(V_Nas2,V_Ks2,0.1.*ones(size(V_Nas2)),'c-','LineWidth',4)


legend('SNIC','HB',...
'Glutamate pulse with Ca-transient',...
'Glutamate pulse, no Ca-transient',...
'Elevated K+ with Ca-transient',...
'Elevated K+, no Ca-transient',...
'Location','northeast');

%focus on one section where the astrocyte trajectory leads through a
%no-oscillation-zone with Ca effects
% load('WB_neuron_ena_ek_loop_Iapp_0165.mat');
freqs = reshape(freq(8,:,:),size(E_Ks,2),size(E_Nas,2));
[x,y] = meshgrid(E_Nas,E_Ks);
x = x(:);
y = y(:);

freqs(freqs==0) = NaN;
figure(20); clf;
scatter(x,y,30,freqs(:),'LineWidth',4); hold on
colorbar
xlabel('E_{Na} (mV)'); ylabel('E_K (mV)');
xlim([min(E_Nas) max(E_Nas)]); ylim([min(E_Ks) max(E_Ks)]);
title(['Frequency (Hz), I_{app} = ' num2str(Iapps(8))]);
set(gca,'FontSize',20);
%highlight "normal" conditions
plot(55.0662, -90.0758, 'go', 'LineWidth',4,'MarkerSize',20);
%highlight elevated K conditions
plot(55.0572, -62.4107, 'ro', 'LineWidth',4,'MarkerSize',20);
%add astrocyte modulatory trajectories
plot(V_Nas,V_Ks,'m-','LineWidth',4)
plot(V_Nas2,V_Ks2,'c-','LineWidth',4)
legend('frequency, Hz',...
     'resting conditions','elevated K^+',...
    'astrocyte modulatory trajectory, with Ca-transient',...
    'astrocyte modulatory trajectory, no Ca-transient');

inds = find(y<= -80);
figure(21); clf;
scatter(x(inds),y(inds),30,freqs(inds),'LineWidth',4); hold on
colorbar
xlabel('E_{Na} (mV)'); ylabel('E_K (mV)');
xlim([min(E_Nas) max(E_Nas)]); ylim([min(E_Ks) -80]);
title(['Frequency (Hz), I_{app} = ' num2str(Iapps(8))]);
set(gca,'FontSize',20);
%highlight "normal" conditions
plot(55.0662, -90.0758, 'go', 'LineWidth',4,'MarkerSize',20);
%highlight elevated K conditions
plot(55.0572, -62.4107, 'ro', 'LineWidth',4,'MarkerSize',20);
%add astrocyte modulatory trajectories
plot(V_Nas,V_Ks,'m-','LineWidth',4)
plot(V_Nas2,V_Ks2,'c-','LineWidth',4)
legend('frequency, Hz',...
     'resting conditions','elevated K^+',...
    'astrocyte modulatory trajectory, with Ca-transient',...
    'astrocyte modulatory trajectory, no Ca-transient');


%focus on one section where the astrocyte trajectory leads through a
%no-oscillation-zone with Glu effects
freqs = reshape(freq(8,:,:),size(E_Ks,2),size(E_Nas,2));
[x,y] = meshgrid(E_Nas,E_Ks);
x = x(:);
y = y(:);

load('Na_K_outs_Glupulses.mat');
K_out = K_outs{1};
Na_out = Na_outs{1};
K_in = 93.2; %mM  - set to make V_K = -90 mV at rest
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest
V_Nas = (R*T/F).*log(Na_out./Na_in).*1e3;
V_Ks = (R*T/F).*log(K_out./K_in).*1e3;

K_out2 = K_outs{2};
Na_out2 = Na_outs{2};
V_Nas2 = (R*T/F).*log(Na_out2./Na_in).*1e3;
V_Ks2 = (R*T/F).*log(K_out2./K_in).*1e3;

freqs(freqs==0) = NaN;
figure(22); clf;
scatter(x,y,30,freqs(:),'LineWidth',4); hold on
colorbar
xlabel('E_{Na} (mV)'); ylabel('E_K (mV)');
xlim([min(E_Nas) max(E_Nas)]); ylim([min(E_Ks) max(E_Ks)]);
title(['Frequency (Hz), I_{app} = ' num2str(Iapps(8))]);
set(gca,'FontSize',20);
%highlight "normal" conditions
plot(55.0662, -90.0758, 'go', 'LineWidth',4,'MarkerSize',20);
%add astrocyte modulatory trajectories
plot(V_Nas,V_Ks,'b-','LineWidth',4)
plot(V_Nas2,V_Ks2,'r-','LineWidth',4)
legend('frequency, Hz',...
     'resting conditions',...
    'astrocyte modulatory trajectory, with Ca-transient',...
    'astrocyte modulatory trajectory, no Ca-transient');

inds = find(y<= -80);
figure(21); clf;
scatter(x(inds),y(inds),30,freqs(inds),'LineWidth',4); hold on
colorbar
xlabel('E_{Na} (mV)'); ylabel('E_K (mV)');
xlim([min(E_Nas) max(E_Nas)]); ylim([min(E_Ks) -80]);
title(['Frequency (Hz), I_{app} = ' num2str(Iapps(8))]);
set(gca,'FontSize',20);
%highlight "normal" conditions
plot(55.0662, -90.0758, 'go', 'LineWidth',4,'MarkerSize',20);
%add astrocyte modulatory trajectories
plot(V_Nas,V_Ks,'b-','LineWidth',4)
plot(V_Nas2,V_Ks2,'r-','LineWidth',4)
legend('frequency, Hz',...
     'resting conditions',...
    'astrocyte modulatory trajectory, with Ca-transient',...
    'astrocyte modulatory trajectory, no Ca-transient');

