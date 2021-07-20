%Test_Full_System.m
%tracks changes in intra- and extra-cellular ion concentrations
%of extracellular ions:
%       Glu, K+, Ca2+, Na+, H+ (Glu_out, K_out, Ca_out, Na_out)
%and intracellular ions:
%       K+, Ca2+, Na+ (K_in, Ca_in, Na_in)
%through EAAT2, Kir4.1, NCX, and NKA
%also keeps track of astrocyte membrane potential Vm

addpath('./src'); close all

global g_EAAT2 g_Kir41 g_NCX rho_NKA R F T

%constants
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

%resting astrocyte concentrations
Na_in = 16.6; %mM - Kirischuk et al. 2012
K_in = 140; %mM - Kirischuk et al. 2012
Ca_in = 73*1e-6; %mM - Kirischuk et al. 2012

%elevated synaptic cleft (extracellular) concentrations
% Glu_out = 0.1; %mM (elevated) - Flanagan et al 2018
Glu_out = 0; %mM (rest) - Flanagan et al 2018
K_out = 12; %12mM (elevated) - textbook/Flanagan et al. 2018
% K_out = 1; %mM (rest)
Na_out = 140; %mM - Kirischuk et al. 2012
Ca_out = 2; %mM - Kirischuk et al. 2012
H_out = (10^(-7.4))*1e3; %mM - Kirischuk et al. 2012

Vm = -80; %mV, resting astrocyte membrane potential Kirischuk et al. 2012

X0 = [Na_in; K_in; Na_out; K_out; Ca_out; Glu_out; Vm];

tspan = [0,1e2];

for jj=1:4
if jj==1
    Ca_signal = @(t) (1/.1).*normpdf(t,15,4).*1e-3; %one peak
elseif jj==2
    Ca_signal = @(t) (1/.1).*normpdf(t,15,4).*1e-3 + (0.4/0.1).*normpdf(t,35,4).*1e-3; %two peak
elseif jj==3
    Ca_signal = @(t) (1/.1).*normpdf(t,15,4).*1e-3 + (0.6/0.1).*normpdf(t,35,6).*1e-3; %sustained peak
else
    Ca_signal = @(t) 0.1e-6 + 0.*t;
end

[t,X] = ode15s(@(t,X) system_eqns_with_ca_signal(t,X,Ca_signal), tspan, X0);

labels = {'[Na^+]_i, mM'; '[K^+]_i, mM';...
    '[Na^+]_e, mM'; '[K^+]_e, mM'; '[Ca^{2+}]_e, mM'; '[Glu^-]_e, mM';...
     'V_m, mV';};
for ii=1:7
    figure(ii);
%     plot(t(1:jj-1),X(1:jj-1,ii),'LineWidth',2) %for euler's method plot

    plot(t,X(:,ii),'LineWidth',2); hold on;
    ylabel(labels(ii)); xlabel('t');
    set(gca,'FontSize',20);
end

figure(8);
Ca_signal = Ca_signal(t);
plot(t,Ca_signal,'LineWidth',2); hold on;
    ylabel('prescribed [Ca^{2+}]_i, mM'); xlabel('t');
    set(gca,'FontSize',20);
    
figure(9);
plot(t,g_NCX.*ncx_current(X(:,3),X(:,1),X(:,5),Ca_signal,X(:,7)),'LineWidth',2);
hold on;
ylabel('I_{NCX}'); xlabel('t'); set(gca,'FontSize',20);
end
% figure(9);
% plot(t,g_Kir41.*kir41_current(X(:,4),X(:,2),X(:,7)),'LineWidth',2); hold on;
% % ylabel('I_{Kir}'); 
% % xlabel('t'); set(gca,'FontSize',20);
% 
% % figure(10); clf;
% plot(t,g_EAAT2.*eaat2_current(X(:,3),X(:,6), X(:,2),X(:,4),H_out,X(:,7)),'LineWidth',2); 
% % ylabel('I_{EAAT2}'); xlabel('t'); set(gca,'FontSize',20);
% 
% % figure(11); clf;
% plot(t,rho_NKA.*nka_current(X(:,1),X(:,4)),'LineWidth',2); 
% % ylabel('I_{NKA}'); xlabel('t'); set(gca,'FontSize',20);
% 
% % figure(12); clf;
% plot(t,g_NCX.*ncx_current(X(:,3),X(:,1),X(:,5),Ca_signal(t),X(:,7)),'LineWidth',2); 
% % ylabel('I_{NCX}'); xlabel('t'); set(gca,'FontSize',20);
% ylabel('I'); xlabel('t'); set(gca,'FontSize',20);
% lgd=legend('Kir4.1 current', 'EAAT2 current', 'NKA current', 'NCX current');
% lgd.Location = 'southeast';

% figure(10); clf;
% VKA = (R*T/F)*log(X(:,4)./X(:,2))*1e3; %mV
% plot(t,VKA,'LineWidth',2); 
% ylabel('E_K'); xlabel('t'); set(gca,'FontSize',20);

