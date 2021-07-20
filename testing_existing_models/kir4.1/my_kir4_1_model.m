% Kir4.1 Channel model
% based on I-V data

clear; close all; clc;
%Kir
%constants
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

%I-V curve fitting params
GKir = 0.03; %nA /mV, membrane conductance - made up to get correct scale
VA1 = -3; %mV
% VA2 = 100; %mV
% VA3 = 50; %mV
VA2 = 90; %mV  %use only 2 VA's to get more consistent graphs across datasets

%Vm for I-V curve
Vm = (-165:5:60); %mV

%external/internal[K+]  - could change dynamically
Ke = 5*1e-3; %M, external [K+]  %
Ki = 130*1e-3; %M, internal [K+]  %from Seifart et al. 2009 (same as Chai 2017)

%functions - note RT/F has units volts
VKA = (R*T/F)*log(Ke/Ki)*1e3; %mV
% I_K = GKir.*(Vm-VKA-VA1)./(1+exp((Vm-VKA-VA2)./VA3)); %units nA
% I_K = @(x,Vm) x(1).*(Vm-VKA-x(2))./(1+exp((Vm-VKA-x(3))./x(4))); %units nA
I_K = @(x,Vm) x(1).*(Vm-VKA-x(2))./(1+exp((Vm-VKA-x(2))./x(3))); %units nA

%-----------Seifart et al., 2009 Data
%fit to Seifart et al Fig1 (A2-inset) Ba2+-sensitive current I-V data
[Vm_data, I_data] = csvimport('kir4_1_data/seifart_etal_2009_fig1inset.csv', 'columns', [1, 2] ,'noHeader', true);
% x0 = [GKir,VA1,VA2,VA3];
x0 = [GKir,VA1,VA2];
Seifart_fit = lsqcurvefit(I_K,x0,Vm_data,I_data);

% figure(1);clf;
plot(Vm, I_K(Seifart_fit,Vm), 'b-','LineWidth',3); hold on
% plot([VKA,VKA],[0,0], 'ks');
% text(VKA-20,0.1,['E_K = ',num2str(VKA,'%4.0f'),'mV'])
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
xlim([-165,60]); ylim([-2,2]);

plot(Vm_data, I_data, 'b^','MarkerSize',10,'LineWidth',2)
% legend('model', 'Seifart et al., 2009')

%---- other figures from Seifart et al., 2009 ------
% %fig1b - normalized to "max" inward current
% % [Vm_data2 I_data2] = csvimport('kir4_1_data/seifart_etal_2009_fig1b_IV.csv', 'columns', [1, 2] ,'noHeader', true);
% plot(Vm_data2, I_data2, '^')
% fig9c - normalized to "max" inward current, but clearly not actually normalized in this way
[Vm_dataS2, I_dataS2] = csvimport('kir4_1_data/seifart_etal_2009_fig9c.csv', 'columns', [1, 2] ,'noHeader', true);
I_dataS2 = (I_dataS2./abs(I_dataS2(1)))*abs(I_data(2));
Seifart_fit2 = lsqcurvefit(I_K,Seifart_fit,Vm_dataS2,I_dataS2);
plot(Vm, I_K(Seifart_fit2,Vm), 'k-','LineWidth',3);
plot(Vm_dataS2, I_dataS2, 'k^')

%-----------Sicca et al., 2016 Data
%external/internal[K+]  - could change dynamically
Ke = 5*1e-3; %M, external [K+]  %
Ki = 120*1e-3; %M, internal [K+]  %from Sicca et al., 2016

%functions - note RT/F has units volts
VKA2 = (R*T/F)*log(Ke/Ki)*1e3; %mV
C = 0.05; %pF
% I_K = @(x,Vm) x(1)*(Vm-VKA2-x(2))./(1+exp((Vm-VKA2-x(3))./x(4))); %units nA
I_K = @(x,Vm) x(1)*(Vm-VKA2-x(2))./(1+exp((Vm-VKA2-x(2))./x(3))); %units nA

%fit to Sicca et al Fig1C, Wild-Type Ba2+-sensitive current I-V data
[Vm_data2, I_data2] = csvimport('kir4_1_data/sicca_etal_2016_fig1c_WT.csv', 'columns', [1, 2] ,'noHeader', true);
I_data2 = I_data2*C; %nA = (pA/pF)*(pF)*(nA/pA)
Sicca_fit = lsqcurvefit(I_K,Seifart_fit,Vm_data2,I_data2);

% figure(2); clf;
plot(Vm, I_K(Sicca_fit,Vm), 'r-','LineWidth',3); hold on
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
xlim([-165,60]); ylim([-2.5,5]);

plot(Vm_data2, I_data2, 'r^','MarkerSize',10,'LineWidth',2)

%with gain-of-function mutation
[Vm_data3, I_data3] = csvimport('kir4_1_data/sicca_etal_2016_fig1c_R18Q.csv', 'columns', [1, 2] ,'noHeader', true);
I_data3 = I_data3*C; %nA = (pA/pF)*(pF)*(nA/pA)
%with simple 2.3x conductance on previous model
plot(Vm, 2.3*I_K(Sicca_fit,Vm), 'm--','LineWidth',3);
plot(Vm_data3, I_data3, 'm^','MarkerSize',10,'LineWidth',2)
%refit model
Sicca_fit2 = lsqcurvefit(I_K,Sicca_fit,Vm_data3,I_data3);
plot(Vm, I_K(Sicca_fit2,Vm), 'm-','LineWidth',3);

%-----------Chai et al, 2017 -  Hippocampus vs. Striatum data--------------
[Vm_data4, I_data4] = csvimport('kir4_1_data/chai_etal_2017_fig1F_hippocampus.csv', 'columns', [1, 2] ,'noHeader', true);
plot(Vm_data4, I_data4, 'g^','MarkerSize',10,'LineWidth',2)
%have to go back and reuse VKA1 for 5mM/130mM extracellular/intracellular K+ concentration
I_K = @(x,Vm) x(1)*(Vm-VKA-x(2))./(1+exp((Vm-VKA-x(2))./x(3))); %units nA
Chai_fit = lsqcurvefit(I_K,Seifart_fit,Vm_data4,I_data4);
plot(Vm, I_K(Chai_fit,Vm), 'g-','LineWidth',3); 

[Vm_data5, I_data5] = csvimport('kir4_1_data/chai_etal_2017_fig1F_striatum.csv', 'columns', [1, 2] ,'noHeader', true);
plot(Vm_data5, I_data5, 'c^','MarkerSize',10,'LineWidth',2)
Chai_fit2 = lsqcurvefit(I_K,Chai_fit,Vm_data5,I_data5);
plot(Vm, 0.8*I_K(Chai_fit,Vm), 'c--','LineWidth',3); 
plot(Vm, I_K(Chai_fit2,Vm), 'c-','LineWidth',3); 
legend({'model 1', 'Seifart et al., 2009','model 2', 'Seifart 2','model 3','WT- Sicca et al., 2016',...
    'model 3 with 2.3x conductance', 'R18Q- Sicca et al., 2016','model 4',...
    'Hippocampus - Chai et al., 2017','model 5','Striatum - Chai et al., 2017',...
    'model 5 with 0.8x conductance', 'model 6'},...
    'Location','northwest')
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

%---------- Plot separate sets of data/curves with about the same Vr--------------
figure(10); clf;
plot(Vm, I_K(Seifart_fit2,Vm), 'k-','LineWidth',3); hold on
plot(Vm_dataS2, I_dataS2, 'k^')
plot(Vm, I_K(Sicca_fit,Vm), 'r-','LineWidth',3); hold on
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
xlim([-165,60]); ylim([-2.5,5]);
plot(Vm_data2, I_data2, 'r^','MarkerSize',10,'LineWidth',2)
plot(Vm, 2.3*I_K(Sicca_fit,Vm), 'm--','LineWidth',3);
plot(Vm_data3, I_data3, 'm^','MarkerSize',10,'LineWidth',2)
plot(Vm, I_K(Sicca_fit2,Vm), 'm-','LineWidth',3);
legend({'model 2', 'Seifart 2','model 3','WT- Sicca et al., 2016',...
    'model 3 with 2.3x conductance', 'R18Q- Sicca et al., 2016','model 4'},...
    'Location','northwest')
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

figure(11); clf;
plot(Vm, I_K(Seifart_fit,Vm), 'b-','LineWidth',3); hold on
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
xlim([-165,60]); ylim([-2,2]);
plot(Vm_data, I_data, 'b^','MarkerSize',10,'LineWidth',2)
plot(Vm_data4, I_data4, 'g^','MarkerSize',10,'LineWidth',2)
plot(Vm, I_K(Chai_fit,Vm), 'g-','LineWidth',3); 
plot(Vm_data5, I_data5, 'c^','MarkerSize',10,'LineWidth',2)
plot(Vm, 0.8*I_K(Chai_fit,Vm), 'c--','LineWidth',3); 
plot(Vm, I_K(Chai_fit2,Vm), 'c-','LineWidth',3); 
legend({'model 1', 'Seifart et al., 2009', ...
    'Hippocampus - Chai et al., 2017','model 5','Striatum - Chai et al., 2017',...
    'model 5 with 0.8x conductance', 'model 6'},...
    'Location','northwest')
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');


%---------- Plot model fit parameters --------------
Gkirs = [Seifart_fit(1), Sicca_fit(1), Sicca_fit2(1), Chai_fit(1), Chai_fit2(1)];
VA1s = [Seifart_fit(2), Sicca_fit(2), Sicca_fit2(2), Chai_fit(2),Chai_fit2(2)];
VA2s = [Seifart_fit(3), Sicca_fit(3), Sicca_fit2(3), Chai_fit(3),Chai_fit2(3)];
figure(2);
subplot(1,3,1); bar(1,Seifart_fit(1),'b'); hold on;
bar(2,Seifart_fit2(1),'k');
bar(3,Sicca_fit(1),'r'); bar(4,Sicca_fit2(1),'m');
bar(5,Chai_fit(1),'g'); bar(6,Chai_fit2(1),'c');
ylabel('nA/mV'); title('G_{Kir}')
set(gca,'XTick',1:6,'xticklabel',{'1','2','3','4','5','6'},'FontSize',20)

subplot(1,3,2); bar(1,Seifart_fit(2),'b'); hold on;
bar(2,Seifart_fit2(2),'k');
bar(3,Sicca_fit(2),'r'); bar(4,Sicca_fit2(2),'m');
bar(5,Chai_fit(2),'g'); bar(6,Chai_fit2(2),'c');
ylabel('mV');title('V_1')
set(gca,'XTick',1:6,'xticklabel',{'1','2','3','4','5','6'},'FontSize',20)

subplot(1,3,3); bar(1,Seifart_fit(3),'b'); hold on;
bar(2,Seifart_fit2(3),'k');
bar(3,Sicca_fit(3),'r'); bar(4,Sicca_fit2(3),'m');
bar(5,Chai_fit(3),'g'); bar(6,Chai_fit2(3),'c');
ylabel('mV'); title('V_2')
set(gca,'XTick',1:6,'xticklabel',{'1','2','3','4','5','6'},'FontSize',20)

%---------- Model at different [K+]_e as in Ransom & Sontheimer, 1995 --------------
Ki = 130*1e-3; %M

%using Seifart_fit
Kes = [2, 3, 5, 10, 20]*1e-3; %M, external [K+] 
%functions - note RT/F has units volts
VKAs = (R*T/F).*log(Kes./Ki).*1e3; %mV
figure(3); clf;
for ii=1:size(Kes,2)
   I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
   plot(Vm, I_K(Seifart_fit,Vm), '-','LineWidth',3); hold on   
end
legend('[K^+]_e = 2 mM','[K^+]_e = 3 mM','[K^+]_e = 5 mM',...
    '[K^+]_e = 10 mM','[K^+]_e = 20 mM','Location','southeast');
xlim([-165,60]); ylim([-3,3]); 
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
title('Model 1 (Seifart et al., 2009 Fit)');

%using Sicca_fit
Kes = [2, 3, 5, 10, 20]*1e-3; %M, external [K+] 
%functions - note RT/F has units volts
VKAs = (R*T/F).*log(Kes./Ki).*1e3; %mV
figure(4); clf;
for ii=1:size(Kes,2)
   I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
   plot(Vm, I_K(Sicca_fit,Vm), '-','LineWidth',3); hold on   
end
legend('[K^+]_e = 2 mM','[K^+]_e = 3 mM','[K^+]_e = 5 mM',...
    '[K^+]_e = 10 mM','[K^+]_e = 20 mM','Location','southeast');
xlim([-165,60]); ylim([-3,3]); 
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
title('Model 2 (Sicca et al., 2016 WT Fit)');

%using Chai_fit
Kes = [2, 3, 5, 10, 20]*1e-3; %M, external [K+] 
%functions - note RT/F has units volts
VKAs = (R*T/F).*log(Kes./Ki).*1e3; %mV
figure(5); clf;
for ii=1:size(Kes,2)
   I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
   plot(Vm, I_K(Chai_fit,Vm), '-','LineWidth',3); hold on   
end
legend('[K^+]_e = 2 mM','[K^+]_e = 3 mM','[K^+]_e = 5 mM',...
    '[K^+]_e = 10 mM','[K^+]_e = 20 mM','Location','southeast');
xlim([-165,60]); ylim([-3,3]); 
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
title('Model 4 (Chai et al., 2017 Hippocampus Fit)');

%--------- Compare directly to Ransom & Sontheimer 1995
%Ransom_Sontheimer data: zero-crossings of whole-cell inward potassium currents,
%not specifically Kir4.1
[Vzero,Izero] = csvimport('kir4_1_data/ransom_sontheimer_fig4b_Zero crossings.csv','columns', [1, 2] ,'noHeader', true);
Vzero = flip(Vzero); Izero=flip(Izero);
Kes = [3,10,20]; %M
Ki = 145; %M as in Ransom & Sontheimer 1995
figure(6); clf;
subplot(1,3,1);
plot(Vzero,Izero,'s','MarkerSize',20); hold on
%using Seifart_fit
%functions - note RT/F has units volts
VKAs = (R*T/F).*log(Kes./Ki).*1e3; %mV
for ii=1:size(Kes,2)
   I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
   plot(Vm, I_K(Seifart_fit,Vm), '-','LineWidth',3); hold on   
end
legend('V_r, R&S 1995','[K^+]_e = 3 mM',...
    '[K^+]_e = 10 mM','[K^+]_e = 20 mM','Location','southeast');
xlim([-165,60]); ylim([-3,3]); 
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
title('Model 1 (Seifart Fit)');

subplot(1,3,2);
plot(Vzero,Izero,'s','MarkerSize',20); hold on
%using Sicca_fit
%functions - note RT/F has units volts
VKAs = (R*T/F).*log(Kes./Ki).*1e3; %mV
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
title('Model 2 (Sicca WT Fit)');

subplot(1,3,3);
plot(Vzero,Izero,'s','MarkerSize',20); hold on
%using Chai_fit
%functions - note RT/F has units volts
VKAs = (R*T/F).*log(Kes./Ki).*1e3; %mV
for ii=1:size(Kes,2)
   I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
   plot(Vm, I_K(Chai_fit,Vm), '-','LineWidth',3); hold on   
end
legend('V_r, R&S 1995','[K^+]_e = 3 mM',...
    '[K^+]_e = 10 mM','[K^+]_e = 20 mM','Location','southeast');
xlim([-165,60]); ylim([-3,3]); 
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
title('Model 4 (Chai HC Fit)');


%--- Plot of Nernst potential E_K as function of [K+]e
Kes = logspace(0,2,10); %mV
Ki = 145; %mV
VKAs = (R*T/F).*log(Kes./Ki).*1e3; %mV
figure(7); clf
semilogx(Kes,VKAs,'s-','LineWidth',2,'MarkerSize',20); hold on;
xlabel('[K^+]_e (mM)'); ylabel('mV');
%--- Plot of Resting membrane potential E_k as function of [K+]e
Vrms1 = zeros(size(VKAs));
Vrms2 = zeros(size(VKAs));
Vrms4 = zeros(size(VKAs));
for ii=1:size(Kes,2)
    I_K = @(x,Vm) x(1)*(Vm-VKAs(ii)-x(2))./(1+exp((Vm-VKAs(ii)-x(2))./x(3))); %units nA
    Vrms1(ii) = fzero(@(Vm) I_K(Seifart_fit,Vm),VKAs(ii)+Seifart_fit(2));
    Vrms2(ii) = fzero(@(Vm) I_K(Sicca_fit,Vm),VKAs(ii)+Sicca_fit(2));
    Vrms4(ii) = fzero(@(Vm) I_K(Chai_fit,Vm),VKAs(ii)+Chai_fit(2));
end
semilogx(Kes,Vrms1,'o--','LineWidth',3,'MarkerSize',10);
semilogx(Kes,Vrms2,'o--','LineWidth',3,'MarkerSize',10);
semilogx(Kes,Vrms4,'o--','LineWidth',3,'MarkerSize',10);
legend('Nernst potential','V_{rest}, model 1 (Seifart fit)',...
    'V_{rest}, model 2 (Sicca WT fit)','V_{rest}, model 4 (Chai HC fit)',...
    'Location','southeast');
set(gca,'FontSize',20);
