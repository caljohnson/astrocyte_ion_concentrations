%wang-buzsaki hippocampal neuron astrocyte effects timingsed current inject
%test to see if faster astrocytic potassium uptake and sodium dynamics have an effect of neural
%excitability

%demonstrates that elevated potassium leads to increased excitability
% and that Ca-transient helps lower excitability to compensate

% load('remote_WB_neuron_ena_ek_loop_v1.mat'); %v1 - EK-95 to -50, ENa50 to 60
% load('remote_WB_neuron_ena_ek_loop_v2.mat'); %v2 - -95 to -80, ENa50 to 60
% load('remote_WB_neuron_ena_ek_loop_v3.mat'); %v3 - EK-95 to -80, ENa50 to 60, Iapp 0.1-0.2
load('remote_WB_neuron_ena_ek_loop_v4.mat'); %v4 - EK-95 to -50, ENa50 to 60, Iapp 0.1-0.2

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

EKmax = -56;
x1 = x(y<EKmax);
y1 = y(y<EKmax);
z1 = z(y<EKmax);

x2 = x(y>=EKmax & z>0.11);
y2 = y(y>=EKmax & z>0.11);
z2 = z(y>=EKmax & z>0.11);


% B = [x1(:) y1(:) ones(size(x1(:)))] \ z1(:);
% xv = linspace(min(x1), max(x1), 100)';
% yv = linspace(min(y1), max(y1), 100)';
% [X1,Y1] = meshgrid(xv, yv);
% Z1 = reshape([X1(:), Y1(:), ones(size(X1(:)))] * B, numel(xv), []);
% mesh(X1, Y1, Z1,'FaceColor','g', 'FaceAlpha', 0.5)
% 
% B2 = [x2(:) y2(:) ones(size(x2(:)))] \ z2(:);
% xv2 = linspace(min(x2), max(x2), 100)';
% yv2 = linspace(min(y2), max(y2), 100)';
% [X2,Y2] = meshgrid(xv2, yv2);
% Z2 = reshape([X2(:), Y2(:), ones(size(X2(:)))] * B2, numel(xv2), []);
% mesh(X2, Y2, Z2,'FaceColor','r', 'FaceAlpha', 0.5)
% 
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
legend('SNIC','HB', 'Location','northeast');