% Kir4.1 Channel model
% based on I-V data from Ransom & Sontheimer 1995 Fig 7b

clear; clc;
%Kir
%constants
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 295.15;%22 C = 295.15 K 

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
VKA = (R*T/F)*log(Ke/Ki)*1e3; %mV

I_K = @(x,Vm) x(1).*(Vm-VKA-x(2))./(1+exp((Vm-VKA-x(2))./x(3))); %units nA

%-----------Ransom & Sontheimer 1995 Data
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

Vm_data = Vm_data_control;
I_data = I_data_control-I_data_ba100;

x0 = [GKir,VA1,VA2];
RS_fit = lsqcurvefit(I_K,x0,Vm_data,I_data);

figure(1);clf;
plot(Vm, I_K(RS_fit,Vm), 'b-','LineWidth',3); hold on
plot([VKA,VKA],[0,0], 'ks');
text(VKA-20,0.1,['E_K = ',num2str(VKA,'%4.0f'),'mV'])
xlabel('membrane potential, mV'); ylabel('current, nA');
set(gca,'FontSize',20);
xlim([-165,60]); ylim([-5,2]);

plot(Vm_data, I_data, 'b^','MarkerSize',10,'LineWidth',2)
legend('model', 'Ransom & Sontheimer 1995');
%axes lines
line([-165,60],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
