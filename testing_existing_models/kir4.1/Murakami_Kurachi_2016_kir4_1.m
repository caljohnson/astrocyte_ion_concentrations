% Murakami & Kurachi, 2016 - Astrocyte Model
% Looking at their astrocyte K+ channel model
% specifically as a model of Kir4.1
% coding up the IK from this paper (equation 29 for I_Kir, p33)

clear; close all; clc;
%Kir
%constants
GKir_syn = 4.4*1e-6; %pA/M^1/2 \mu.m^2 mV, max membrane conductance
GKir_main = 4.4*1e-7; %pA/M^1/2 \mu.m^2 mV, max membrane conductance
GKir_peri = 4.4*1e-5; %pA/M^1/2 \mu.m^2 mV, max membrane conductance
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

%fitted params
A_syn = 6750;   %\mu m^2, membrane area
A_main = 37464.6;  %\mu m^2, membrane area
A_peri = 785.4; %\mu m^2, membrane area
VA2 = 91.7; %mV
VA3 = 25.6; %mV

%Vm for I-V curve
Vm = (-140:5:60); %mV

%external/internal[K+]  - could change dynamically
% Ke = 2.5*1e-3; %M, external [K+]  %
% Ki = 63*1e-3; %M, internal [K+]  %picked so that VKA = -86 mV
%initial, in paper
Ke = 3e-3; %M, external [K+]
Ki = 130*1e-3; %M, internal [K+] -> gives VKA = -100.6 mV, not correct reversal?
%model technically has different Ke,Ki in each compartment, but not clear
%where to find these values - vary dynamically?

%functions
VKA = (R*T/F)*log(Ke/Ki)*1e3; %mV
I_K_syn = GKir_syn*A_syn.*(Vm-VKA).*sqrt(Ke)./(1+exp((Vm-VKA-VA2)./VA3)); %units pA
I_K_main = GKir_main*A_main.*(Vm-VKA).*sqrt(Ke)./(1+exp((Vm-VKA-VA2)./VA3)); %units pA
I_K_peri = GKir_peri*A_peri.*(Vm-VKA).*sqrt(Ke)./(1+exp((Vm-VKA-VA2)./VA3)); %units pA
I_K = I_K_syn+I_K_main+I_K_peri;

plot(Vm,I_K_syn); hold on
plot(Vm,I_K_main);
plot(Vm,I_K_peri);
plot(Vm, I_K, '-o'); hold on
xlabel('membrane potential, mV'); ylabel('current, pA');
set(gca,'FontSize',20);
xlim([-140,60]); %ylim([-200,1400]);
line([-140,60],[0,0],'LineStyle','--', 'Color', 'k');
line([0,0],[-.3,0.3],'LineStyle','--', 'Color', 'k');
return

%external/internal[K+]  - could change dynamically
Ke = 3*1e-3; %M, external [K+]  %bath solution, from Seifart et al., 2009
Ki = 145*1e-3; %M, internal [K+]  %pipette solution, from Seifart et al., 2009
%functions
%note RT/F has units volts
VKA = (R*T/F)*log(Ke/Ki)*1e3; %mV
I_K = GKir*A.*(Vm-VKA).*sqrt(Ke)./(1+exp((Vm-VKA-VA2)./VA3)); %units pA
%plot
plot(Vm, I_K, '-o'); hold on


%external/internal[K+]  - could change dynamically
Ke = 5*1e-3; %M, external [K+]  %bath solution, from Seifart et al., 2009
Ki = 135*1e-3; %M, internal [K+]  %pipette solution, from Seifart et al., 2009
%functions
VKA = (R*T/F)*log(Ke/Ki)*1e3; %mV
I_K = GKir*A.*(Vm-VKA).*sqrt(Ke)./(1+exp((Vm-VKA-VA2)./VA3)); %units pA
%plot
plot(Vm, I_K, '-o'); hold on

%external/internal[K+]  - could change dynamically
Ke = 10*1e-3; %M, external [K+]  %bath solution, from Seifart et al., 2009
Ki = 135*1e-3; %M, internal [K+]  %pipette solution, from Seifart et al., 2009
%functions
VKA = (R*T/F)*log(Ke/Ki)*1e3; %mV
I_K = GKir*A.*(Vm-VKA).*sqrt(Ke)./(1+exp((Vm-VKA-VA2)./VA3)); %units pA
%plot
plot(Vm, I_K, '-o'); hold on



legend('Ke/Ki = 2.5/135', 'Ke/Ki = 3/145', 'Ke/Ki = 5/135', 'Ke/Ki = 10/135');
line([-140,60],[0,0],'LineStyle','--', 'Color', 'k');
line([0,0],[-4e-4,2e-4],'LineStyle','--', 'Color', 'k');
