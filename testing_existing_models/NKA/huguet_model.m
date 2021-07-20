%NKA model - Huguet et al. 2016

p=1;
kK = 2; %mM, half activation [K+]e concentration for the astrocytic NKA
kNa = 7.7; %mM, half activation [Na+]i concentration for the astrocytic NKA
I_NKA = @(K_out, Na_in) p.*(K_out./(kK + K_out)).^2.*(Na_in./(kNa + Na_in)).^3;

%constants
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 298.15; %K, absolute temperature (25 C in Kirischuk et al. 2012)

%cellular concentrations - Kirischuk et al. 2012
%resting concentrations
Na_in = 16.6; %mM
Ca_in = 73*1e-6; %mM
Na_out = 140; %mM
Ca_out = 2; %mM
K_out = 5;%mM
K_in = 140;%mM

K_outs = linspace(0,100);
figure(1);clf
plot(K_outs,I_NKA(K_outs,Na_in),'b-','LineWidth',4); hold on;
% line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
% line([0,0],[-2,2],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
xlabel('[K^+]_e (mM)'); ylabel('current');set(gca,'FontSize',20);

Na_ins = linspace(0,100);
figure(2);clf
plot(Na_ins,I_NKA(K_out,Na_ins),'b-','LineWidth',4); hold on
xlabel('[Na^+]_i (mM)'); ylabel('current');set(gca,'FontSize',20);

