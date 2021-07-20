% Huguet et al., 2014 - Astrocyte Model
% Looking at their astrocyte K+ channel model
% specifically as a model of Kir4.1
% coding up the IK from this paper (equation for I_K,A p 453)

%constants
PK = 4.810*1e-6; %cm/s, K+ membrane permeability constant
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

%parameters - could change dynamically
Ke = 5*1e-3; %M, external [K+]  %bath solution, from Seifart et al., 2009
Ki = 130*1e-3; %M, internal [K+]  %pipette solution, from Seifart et al., 2009

%Vm for I-V curve
Vm = (-140:5:60).*1e-3; %V, wheras I-V curves in other papers are in mV

%functions
phi = Vm.*F./(R*T); %dimensionless since F/RT has units 1/V
I_K = PK*F*phi.*(Ke*exp(-phi) - Ki)./(exp(-phi)-1); %units C/cm^2 s or amp/cm^2

plot(Vm*10^3, I_K*1e4, '-o');
line([-140,60],[0,0],'LineStyle','--', 'Color', 'k');
line([0,0],[-200,1400],'LineStyle','--', 'Color', 'k');
xlabel('membrane potential, mV'); ylabel('current density, A/m^2');
set(gca,'FontSize',20); xlim([-140,60]); ylim([-200,1400]);
%units mV, amp/m^2

% %divide by membrane capacitance to get I = V/s
% Cma = 1e-6*(10^4); %1 muF/cm^2, astrocyte membrane capacitance per area
% plot(Vm*10^3, I_K*1e4./Cma, '-o');
% line([-140,60],[0,0],'LineStyle','--', 'Color', 'k');
% % line([0,0],[-200,1400],'LineStyle','--', 'Color', 'k');
% xlabel('membrane potential, mV'); ylabel('volts per second (A/F)');
% set(gca,'FontSize',20); xlim([-140,60]); %ylim([-200,1400]);
% 

