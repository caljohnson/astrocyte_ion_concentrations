%wang-buzsaki hippocampal neuron astrocyte effects timingsed current inject
%test to see if faster astrocytic potassium uptake and sodium dynamics have an effect of neural
%excitability

load('remote_WB_neuron_gna_gk_loop_v1.mat'); %v1 

[x,y] = meshgrid(g_Nas,g_Ks);
x = x(:);
y = y(:);

for ii=1:size(Iapps,2)
         
    freqs = reshape(freq(ii,:,:),size(g_Ks,2),size(g_Nas,2));
    delays = reshape(delay(ii,:,:),size(g_Ks,2),size(g_Nas,2));
   
    freqs(freqs==0) = NaN;
    figure(ii);
%     surf(E_Nas,E_Ks,freqs); view(2); shading flat; 
    scatter(x,y,30,freqs(:),'LineWidth',4); hold on
    colorbar
xlabel('g_{Na} (mS/cm^2)'); ylabel('g_K (mS/cm^2)');
xlim([min(g_Nas) max(g_Nas)]); ylim([min(g_Ks) max(g_Ks)]);
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

    xlabel('g_{Na} (mS/cm^2)'); ylabel('g_K (mS/cm^2)'); zlabel('I_{app} (mA)');
    zlim([0.1,0.3]);
    
    x = [x; x_bif{ii}];
    y = [y; y_bif{ii}];
    z = [z; I_bif{ii}];

end

gKmax = 4;
gNamin = 40;
x1 = x(y<=gKmax & x>gNamin);
y1 = y(y<=gKmax & x>gNamin);
z1 = z(y<=gKmax & x>gNamin);

x2 = x(y>gKmax);
y2 = y(y>gKmax);
z2 = z(y>gKmax);


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

return
 