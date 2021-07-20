% Sibille et al., 2015 - Astrocyte Model
% Looking at their astrocyte K+ channel model
% specifically as a model of Kir4.1
% coding up the IK from this paper (equation 22 for I_Kir, p16)

clear; clc; close all;
%constants
GKir = 60*1e-12; %S, 60pS-> Kir membrane conductance (Siemen = amp/V)
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

%fitted params
VA1 = -14.83*1e-3;%V
VA2 = 34*1e-3;    %V
VA3 = 19.2*1e-3; %V

%Vm for I-V curve
Vm = (-140:5:60).*1e-3; %V, wheras I-V curves in other papers are in mV

%external/internal[K+]  - could change dynamically
Ke = 2.5*1e-3; %M, external [K+]  %bath solution, from Seifart et al., 2009
Ki = 135*1e-3; %M, internal [K+]  %pipette solution, from Seifart et al., 2009

%functions
VKA = (R*T/F)*log(Ke/Ki); %V
I_K = GKir.*(Vm-VKA-VA1).*sqrt(Ke)./(1+exp((Vm-VKA-VA2)./VA3)); %units amp*M^(1/2)

plot(Vm*10^3, I_K*1e9, '-o'); hold on
xlabel('membrane potential, mV'); ylabel('current times Molarity^{1/2} , nA*M^{1/2}');
set(gca,'FontSize',20);
xlim([-140,60]); %ylim([-200,1400]);

%external/internal[K+]  - could change dynamically
Ke = 3*1e-3; %M, external [K+]  %bath solution, from Seifart et al., 2009
Ki = 145*1e-3; %M, internal [K+]  %pipette solution, from Seifart et al., 2009
%functions
VKA = (R*T/F)*log(Ke/Ki); %V
I_K = GKir.*(Vm-VKA-VA1).*sqrt(Ke)./(1+exp((Vm-VKA-VA2)./VA3)); %units amp*M^(1/2)
plot(Vm*10^3, I_K*1e9, '-o');

%external/internal[K+]  - could change dynamically
Ke = 5*1e-3; %M, external [K+]  %bath solution, from Seifart et al., 2009
Ki = 135*1e-3; %M, internal [K+]  %pipette solution, from Seifart et al., 2009
%functions
VKA = (R*T/F)*log(Ke/Ki); %V
I_K = GKir.*(Vm-VKA-VA1).*sqrt(Ke)./(1+exp((Vm-VKA-VA2)./VA3)); %units amp*M^(1/2)
plot(Vm*10^3, I_K*1e9, '-o');

%external/internal[K+]  - could change dynamically
Ke = 10*1e-3; %M, external [K+]  %bath solution, from Seifart et al., 2009
Ki = 135*1e-3; %M, internal [K+]  %pipette solution, from Seifart et al., 2009
%functions
VKA = (R*T/F)*log(Ke/Ki); %V
I_K = GKir.*(Vm-VKA-VA1).*sqrt(Ke)./(1+exp((Vm-VKA-VA2)./VA3)); %units amp*M^(1/2)
plot(Vm*10^3, I_K*1e9, '-o');
legend('Ke/Ki = 2.5/135', 'Ke/Ki = 3/145', 'Ke/Ki = 5/135', 'Ke/Ki = 10/135');

line([-140,60],[0,0],'LineStyle','--', 'Color', 'k');
% line([0,0],[-4e-4,2e-4],'LineStyle','--', 'Color', 'k');
