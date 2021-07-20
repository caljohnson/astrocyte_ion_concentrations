% Kir4.1 Channel model
% based on I-V data

clear; clc;
%Kir
%constants
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 295.15;%22 C = 295.15 K 
%T = 310; %K, absolute temperature

%I-V curve fitting params
GKir = 0.03; %nA /mV, membrane conductance - made up to get correct scale
VA1 = -3; %mV
% VA2 = 100; %mV
% VA3 = 50; %mV
VA2 = 90; %mV  %use only 2 VA's to get more consistent graphs across datasets

%Vm for I-V curve
Vm = (-165:5:60); %mV

Inorm = 1.6; %approximate current at 0mV from Chai 2017 Hippocampus data

%external/internal[K+]  - could change dynamically
Ke = 5*1e-3; %M, external [K+]  %from Seifart et al. 2009
Ki = 130*1e-3; %M, internal [K+]  %from Seifart et al. 2009

%functions - note RT/F has units volts
VKA1 = (R*T/F)*log(Ke/Ki)*1e3; %mV
I_K1 = @(x,Vm) x(1).*(Vm-VKA1-x(2))./(1+exp((Vm-VKA1-x(2))./x(3))); %units nA


%---- data from Seifart et al., 2009 ------
% %fig9C - normalized to "max" inward current
% fig9c - normalized to "max" inward current, but clearly not actually normalized in this way
[Vm_data1, I_data2] = csvimport('kir4_1_data/seifart_etal_2009_fig9c.csv', 'columns', [1, 2] ,'noHeader', true);
I_data2 = (I_data2./abs(I_data2(end-3)))*Inorm; %normalized to match Chai Hippocampus at 0mV
x0 = [GKir,VA1,VA2];
Seifart_fit = lsqcurvefit(I_K1,x0,Vm_data1,I_data2);
fig=figure(1);clf;
plot(Vm, I_K1(Seifart_fit,Vm), 'k-','LineWidth',3); hold on
plot(Vm_data1, I_data2, 'k^','MarkerSize',10,'LineWidth',2)
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
xlim([-165,60]); ylim([-2,2]);

%-----------Ransom & Sontheimer 1995 Data
%external/internal[K+]  - could change dynamically
Ke = 3*1e-3; %M, external [K+]  %from R&S 1995
Ki = 145*1e-3; %M, internal [K+]  %from R&S 1995

%functions - note RT/F has units volts
VKA4 = (R*T/F)*log(Ke/Ki)*1e3; %mV
I_K1 = @(x,Vm) x(1).*(Vm-VKA4-x(2))./(1+exp((Vm-VKA4-x(2))./x(3))); %units nA

%fit to Fig7B Ba2+-sensitive current I-V data
[Vm_data_control, I_data_control] = ...
    csvimport('kir4_1_data/ransom_sontheimer_1995_fig7b_control.csv',...
    'columns', [1, 2] ,'noHeader', true);
%get rid of extra entry
Vm_data_control = [Vm_data_control(1:5);Vm_data_control(7:end);];
I_data_control = [I_data_control(1:5);I_data_control(7:end);];

[Vm_data_ba100, I_data_ba100] = ...
    csvimport('kir4_1_data/ransom_sontheimer_1995_fig7b_ba100.csv',...
    'columns', [1, 2] ,'noHeader', true);

RS_fit= lsqcurvefit(I_K1,Seifart_fit,Vm_data_control,I_data_control-I_data_ba100);
plot(Vm, I_K1(RS_fit,Vm), '-','LineWidth',3,...
    'Color',[1.0000    0.4118    0.1608]); hold on
plot(Vm_data_control, I_data_control-I_data_ba100, '^','MarkerSize',10,'LineWidth',2,...
    'Color',[1.0000    0.4118    0.1608]); 

%-----------Sicca et al., 2016 Data
%external/internal[K+]  - could change dynamically
Ke = 5*1e-3; %M, external [K+]  %
Ki = 120*1e-3; %M, internal [K+]  %from Sicca et al., 2016

%functions - note RT/F has units volts
VKA2 = (R*T/F)*log(Ke/Ki)*1e3; %mV
C = 0.05; %pF
% I_K = @(x,Vm) x(1)*(Vm-VKA2-x(2))./(1+exp((Vm-VKA2-x(3))./x(4))); %units nA
I_K2 = @(x,Vm) x(1)*(Vm-VKA2-x(2))./(1+exp((Vm-VKA2-x(2))./x(3))); %units nA

%fit to Sicca et al Fig1C, Wild-Type Ba2+-sensitive current I-V data
[Vm_data2, I_data2] = csvimport('kir4_1_data/sicca_etal_2016_fig1c_WT.csv', 'columns', [1, 2] ,'noHeader', true);
I2max = abs(I_data2(end-3));
I_data2 = I_data2./abs(I_data2(end-3))*Inorm; %nA = (pA/pF)*(pF)*(nA/pA)
Sicca_fit = lsqcurvefit(I_K2,Seifart_fit,Vm_data2,I_data2);

plot(Vm, I_K2(Sicca_fit,Vm), 'r-','LineWidth',3); hold on
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
xlim([-165,60]); ylim([-2.5,5]);

plot(Vm_data2, I_data2, 'r^','MarkerSize',10,'LineWidth',2)

%with gain-of-function mutation
[Vm_data3, I_data3] = csvimport('kir4_1_data/sicca_etal_2016_fig1c_R18Q.csv', 'columns', [1, 2] ,'noHeader', true);
I_data3 = I_data3/I2max*Inorm; %nA = (pA/pF)*(pF)*(nA/pA)
%with simple 2.3x conductance on previous model
plot(Vm, 2.3*I_K2(Sicca_fit,Vm), 'm--','LineWidth',3);
plot(Vm_data3, I_data3, 'm^','MarkerSize',10,'LineWidth',2)
%refit model
Sicca_fit2 = lsqcurvefit(I_K2,Sicca_fit,Vm_data3,I_data3);
plot(Vm, I_K2(Sicca_fit2,Vm), 'm-','LineWidth',3);

%-----------Chai et al, 2017 -  Hippocampus vs. Striatum data--------------
[Vm_data4, I_data4] = csvimport('kir4_1_data/chai_etal_2017_fig1F_hippocampus.csv', 'columns', [1, 2] ,'noHeader', true);
plot(Vm_data4, I_data4, 'c^','MarkerSize',10,'LineWidth',2)
%external/internal[K+]  - could change dynamically
Ke = 4.5*1e-3; %M, external [K+]  %from Chai et al. 2017, given as a concentration on "brain slice artifical CSF"
Ki = 130*1e-3; %M, internal [K+]  %from Chai et al. 2017, given.
%functions - note RT/F has units volts
VKA3 = (R*T/F)*log(Ke/Ki)*1e3; %mV

I_K3 = @(x,Vm) x(1)*(Vm-VKA3-x(2))./(1+exp((Vm-VKA3-x(2))./x(3))); %units nA
Chai_fit = lsqcurvefit(I_K3,Seifart_fit,Vm_data4,I_data4);
plot(Vm, I_K3(Chai_fit,Vm), 'c-','LineWidth',3); 

[Vm_data5, I_data5] = csvimport('kir4_1_data/chai_etal_2017_fig1F_striatum.csv', 'columns', [1, 2] ,'noHeader', true);
plot(Vm_data5, I_data5, 'b^','MarkerSize',10,'LineWidth',2)
Chai_fit2 = lsqcurvefit(I_K3,Chai_fit,Vm_data5,I_data5);
plot(Vm, 0.8*I_K3(Chai_fit,Vm), 'b--','LineWidth',3); 
plot(Vm, I_K3(Chai_fit2,Vm), 'b-','LineWidth',3); 
% legend({'model 1', 'Seifart et al., 2009','model 2', 'Seifart 2','model 3','WT- Sicca et al., 2016',...
%     'model 3 with 2.3x conductance', 'R18Q- Sicca et al., 2016','model 4',...
%     'Hippocampus - Chai et al., 2017','model 5','Striatum - Chai et al., 2017',...
%     'model 5 with 0.8x conductance', 'model 6'},...
%     'Location','northwest')
legend({'model 1', 'Seifart et al., 2009','model 2','Ransom & Sontheimer 1995',...
    'model 3', 'WT- Sicca et al., 2016',...
    'model 3 with 2.3x conductance', 'R18Q- Sicca et al., 2016','model 4',...
    'Hippocampus - Chai et al., 2017','model 5','Striatum - Chai et al., 2017',...
    'model 5 with 0.8x conductance', 'model 6'},...
    'Location','northwest')
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
% saveas(fig,'figs/paper_figs/fig_model_fits_all.png');

%---------- Plot separate sets of data/curves with about the same Vr--------------
fig=figure(2); clf;
% plot(Vm, I_K(Seifart_fit2,Vm), 'k-','LineWidth',3); hold on
% plot(Vm_dataS2, I_dataS2, 'k^','MarkerSize',10,'LineWidth',2)
plot(Vm, I_K2(Sicca_fit,Vm), 'r-','LineWidth',3); hold on 
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
xlim([-165,60]); ylim([-2.5,5]);
plot(Vm_data2, I_data2, 'r^','MarkerSize',10,'LineWidth',2)
plot(Vm, 2.3*I_K2(Sicca_fit,Vm), 'm--','LineWidth',3);
plot(Vm_data3, I_data3, 'm^','MarkerSize',10,'LineWidth',2)
plot(Vm, I_K2(Sicca_fit2,Vm), 'm-','LineWidth',3);
legend({'model 2','WT- Sicca et al., 2016',...
    'model 2 with 2.3x conductance', 'R18Q- Sicca et al., 2016','model 3'},...
    'Location','northwest')
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
% saveas(fig,'figs/paper_figs/fig_model_fits_sicca.png');

fig=figure(3); clf;
plot(Vm, I_K3(Chai_fit,Vm), 'c-','LineWidth',3); hold on;
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
xlim([-165,60]); ylim([-2,2]);
plot(Vm_data4, I_data4, 'c^','MarkerSize',10,'LineWidth',2)
plot(Vm_data5, I_data5, 'b^','MarkerSize',10,'LineWidth',2)
plot(Vm, 0.8*I_K3(Chai_fit,Vm), 'b--','LineWidth',3); 
plot(Vm, I_K3(Chai_fit2,Vm), 'b-','LineWidth',3); 
legend({'model 4', 'Hippocampus - Chai et al., 2017','Striatum - Chai et al., 2017',...
    'model 4 with 0.8x conductance', 'model 5'},...
    'Location','northwest')
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
% saveas(fig,'figs/paper_figs/fig_model_fits_chai.png');

%---------- Plot model fit parameters --------------
Gkirs = [Seifart_fit(1), RS_fit(1), Sicca_fit(1), Sicca_fit2(1), Chai_fit(1), Chai_fit2(1)];
VA1s = [Seifart_fit(2), RS_fit(2), Sicca_fit(2), Sicca_fit2(2), Chai_fit(2),Chai_fit2(2)];
VA2s = [Seifart_fit(3), RS_fit(3), Sicca_fit(3), Sicca_fit2(3), Chai_fit(3),Chai_fit2(3)];
fig=figure(4); clf;
subplot(1,3,1); bar(1,Seifart_fit(1),'k'); hold on;
bar(2,RS_fit(1),'FaceColor',[1.0000    0.4118    0.1608]);
bar(3,Sicca_fit(1),'r'); bar(4,Sicca_fit2(1),'m');
bar(5,Chai_fit(1),'c'); bar(6,Chai_fit2(1),'b');
ylabel('nA/mV'); title('g'); xlabel('model number');
set(gca,'XTick',1:6,'xticklabel',{'1','2','3','4','5','6'},'FontSize',20)

subplot(1,3,2); bar(1,Seifart_fit(2),'k'); hold on;
bar(2,RS_fit(2),'FaceColor',[1.0000    0.4118    0.1608]);
bar(3,Sicca_fit(2),'r'); bar(4,Sicca_fit2(2),'m');
bar(5,Chai_fit(2),'c'); bar(6,Chai_fit2(2),'b');
ylabel('mV');title('V_1'); xlabel('model number');
set(gca,'XTick',1:6,'xticklabel',{'1','2','3','4','5','6'},'FontSize',20)

subplot(1,3,3); bar(1,Seifart_fit(3),'k'); hold on;
bar(2,RS_fit(3),'FaceColor',[1.0000    0.4118    0.1608]);
bar(3,Sicca_fit(3),'r'); bar(4,Sicca_fit2(3),'m');
bar(5,Chai_fit(3),'c'); bar(6,Chai_fit2(3),'b');
ylabel('mV'); title('V_2'); xlabel('model number');
set(gca,'XTick',1:6,'xticklabel',{'1','2','3','4','5','6'},'FontSize',20)
% saveas(fig,'figs/paper_figs/fig_model_params.png');


%---------- Model at different [K+]_e as in Ransom & Sontheimer, 1995 --------------

%using Seifart_fit
Ki = 130*1e-3; %M
Kes = [2, 3, 5, 10, 20]*1e-3; %M, external [K+] 
%functions - note RT/F has units volts
VKAs = (R*T/F).*log(Kes./Ki).*1e3; %mV
fig=figure(5); clf;
subplot(2,2,1);
for ii=1:size(Kes,2)
   I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
   plot(Vm, I_K(Seifart_fit,Vm), '-','LineWidth',3); hold on   
end
lgd = legend('[K^+]_e = 2 mM','[K^+]_e = 3 mM','[K^+]_e = 5 mM',...
    '[K^+]_e = 10 mM','[K^+]_e = 20 mM','Location','northwest');
xlim([-165,60]); ylim([-3,3]); 
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
lgd.FontSize = 15;
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
title('Model 1');

%using RS_fit
Ki = 145*1e-3; %M
Kes = [2, 3, 5, 10, 20]*1e-3; %M, external [K+] 
%functions - note RT/F has units volts
VKAs = (R*T/F).*log(Kes./Ki).*1e3; %mV
subplot(2,2,2);
for ii=1:size(Kes,2)
   I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
   plot(Vm, I_K(RS_fit,Vm), '-','LineWidth',3); hold on   
end
lgd=legend('[K^+]_e = 2 mM','[K^+]_e = 3 mM','[K^+]_e = 5 mM',...
    '[K^+]_e = 10 mM','[K^+]_e = 20 mM','Location','northwest');
xlim([-165,60]); ylim([-3,3]); 
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
lgd.FontSize = 15;
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
title('Model 2');

%using Sicca_fit
Ki = 120*1e-3; %M
Kes = [2, 3, 5, 10, 20]*1e-3; %M, external [K+] 
%functions - note RT/F has units volts
VKAs = (R*T/F).*log(Kes./Ki).*1e3; %mV
subplot(2,2,3);
for ii=1:size(Kes,2)
   I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
   plot(Vm, I_K(Sicca_fit,Vm), '-','LineWidth',3); hold on   
end
lgd= legend('[K^+]_e = 2 mM','[K^+]_e = 3 mM','[K^+]_e = 5 mM',...
    '[K^+]_e = 10 mM','[K^+]_e = 20 mM','Location','northwest');
xlim([-165,60]); ylim([-3,3]); 
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
lgd.FontSize = 15;
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
title('Model 3');

%using Chai_fit
Ki = 130*1e-3; %M
Kes = [2, 3, 5, 10, 20]*1e-3; %M, external [K+] 
%functions - note RT/F has units volts
VKAs = (R*T/F).*log(Kes./Ki).*1e3; %mV
subplot(2,2,4);
for ii=1:size(Kes,2)
   I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
   plot(Vm, I_K(Chai_fit,Vm), '-','LineWidth',3); hold on   
end
lgd = legend('[K^+]_e = 2 mM','[K^+]_e = 3 mM','[K^+]_e = 5 mM',...
    '[K^+]_e = 10 mM','[K^+]_e = 20 mM','Location','northwest');
xlim([-165,60]); ylim([-3,3]); 
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
lgd.FontSize = 15;
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
title('Model 5');

% saveas(fig,'figs/paper_figs/fig_model_Ke_changes.png');

%--------- Compare directly to Ransom & Sontheimer 1995
%Ransom_Sontheimer data: zero-crossings of whole-cell inward potassium currents,
%not specifically Kir4.1
[Vzero,~] = csvimport('kir4_1_data/ransom_sontheimer_fig4b_Zero crossings.csv','columns', [1, 2] ,'noHeader', true);
Vzero = flip(Vzero); Izero=zeros(3,1);
Kes = [3,10,20]; %M
Ki = 145; %M as in Ransom & Sontheimer 1995
fig=figure(6); clf;
subplot(2,2,1);
plot(Vzero,Izero,'bo','MarkerSize',10,'LineWidth',4,'MarkerFaceColor','b'); hold on
%using Seifart_fit
%functions - note RT/F has units volts
VKAs = (R*T/F).*log(Kes./Ki).*1e3; %mV
for ii=1:size(Kes,2)
   I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
   plot(Vm, I_K(Seifart_fit,Vm), '-','LineWidth',3); hold on   
end
legend('V_{rest}, R&S 1995','[K^+]_e = 3 mM',...
    '[K^+]_e = 10 mM','[K^+]_e = 20 mM','Location','southeast');
xlim([-165,60]); ylim([-3,3]); 
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
title('Model 1');

subplot(2,2,2);
plot(Vzero,Izero,'bo','MarkerSize',10,'LineWidth',4,'MarkerFaceColor','b'); hold on
%using RS fit
%functions - note RT/F has units volts
VKAs = (R*T/F).*log(Kes./Ki).*1e3; %mV
for ii=1:size(Kes,2)
   I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
   plot(Vm, I_K(RS_fit,Vm), '-','LineWidth',3); hold on   
end
legend('V_{rest}, R&S 1995','[K^+]_e = 3 mM',...
    '[K^+]_e = 10 mM','[K^+]_e = 20 mM','Location','southeast');
xlim([-165,60]); ylim([-3,3]); 
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
title('Model 2');

subplot(2,2,3);
plot(Vzero,Izero,'bo','MarkerSize',10,'LineWidth',4,'MarkerFaceColor','b'); hold on
%using Sicca_fit
for ii=1:size(Kes,2)
   I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
   plot(Vm, I_K(Sicca_fit,Vm), '-','LineWidth',3); hold on   
end
legend('V_r, R&S 1995','[K^+]_e = 3 mM',...
    '[K^+]_e = 10 mM','[K^+]_e = 20 mM','Location','southeast');
xlim([-165,60]); ylim([-3,3]); 
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
title('Model 3');

subplot(2,2,4);
plot(Vzero,Izero,'bo','MarkerSize',10,'LineWidth',4,'MarkerFaceColor','b'); hold on
%using Chai_fit
for ii=1:size(Kes,2)
   I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
   plot(Vm, I_K(Chai_fit,Vm), '-','LineWidth',3); hold on   
end
legend('V_{rest}, R&S 1995','[K^+]_e = 3 mM',...
    '[K^+]_e = 10 mM','[K^+]_e = 20 mM','Location','southeast');
xlim([-165,60]); ylim([-3,3]); 
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
title('Model 5');
% saveas(fig,'figs/paper_figs/fig_model_Ke_chances_with_Vr.png');

%--- Plot of Nernst potential E_K as function of [K+]e
Kes = logspace(0,2,10); %mV
Ki = 145; %mV
VKAs = (R*T/F).*log(Kes./Ki).*1e3; %mV
fig=figure(7); clf
semilogx(Kes,VKAs,'s-','LineWidth',4,'MarkerSize',20,'Color',[0.9294    0.6941    0.1255]); hold on;
xlabel('[K^+]_e (mM)'); ylabel('mV');
%--- Plot of Resting membrane potential E_k as function of [K+]e
Vrms1 = zeros(size(VKAs));
Vrms2 = zeros(size(VKAs));
Vrms3 = zeros(size(VKAs));
Vrms4 = zeros(size(VKAs));
for ii=1:size(Kes,2)
    I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
    Vrms1(ii) = fzero(@(Vm) I_K(Seifart_fit,Vm),VKAs(ii)+Seifart_fit(2));
    Vrms2(ii) = fzero(@(Vm) I_K(RS_fit,Vm),VKAs(ii)+RS_fit(2));
    Vrms3(ii) = fzero(@(Vm) I_K(Sicca_fit,Vm),VKAs(ii)+Sicca_fit(2));
    Vrms4(ii) = fzero(@(Vm) I_K(Chai_fit,Vm),VKAs(ii)+Chai_fit(2));
end
semilogx(Kes,Vrms1,'ko--','LineWidth',3,'MarkerSize',10);
semilogx(Kes,Vrms2,'o--','LineWidth',3,'MarkerSize',10,...
    'Color',[1.0000    0.4118    0.1608]);
semilogx(Kes,Vrms3,'ro--','LineWidth',3,'MarkerSize',10);
semilogx(Kes,Vrms4,'bo--','LineWidth',3,'MarkerSize',10);
legend('Nernst potential','V_{rest}, model 1',...
    'V_{rest}, model 2','V_{rest}, model 3','V_{rest}, model 5',...
    'Location','southeast');
set(gca,'FontSize',20);
% saveas(fig,'figs/paper_figs/fig_model_Vr_changes.png');
