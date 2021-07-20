% NCX - Sodium-Calcium exchanger model
% from Gall & Susa 1999
% still need to fit to astrocyte data


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

%functions - note RT/F has units volts
%reversal potential for NCX
V_NaCa = (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out/Ca_in))*1e3; %mV
V_NaCa;

%cellular concentrations - Kirischuk et al. 2012
%stimulated concentrations
Na_in = 25; %mM
Ca_in = 4*1e-3; %mM
Na_out = 140; %mM
Ca_out = 2; %mM

%functions - note RT/F has units volts
%reversal potential for NCX
V_NaCa = (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out/Ca_in))*1e3; %mV
V_NaCa;

%hill functions parameters
nH = 1;
K_H = 1.5*1e-3; %1.5 microM
x = 5; %pS
V_NaCa = @(Ca_in) (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out/Ca_in))*1e3; %mV
I_NaCa = @(Vm,Ca_in) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa(Ca_in));

%Vm for I-V curve
Vm = (-100:5:50); %mV

Ca_ins = [600e-6, 500e-6, 400e-6, 300e-6, 200e-6, 100e-6];
%[73e-6, 73e-5,73e-4,4e-4,4e-3,1.5e-3,4e-2];
fig=figure(1);clf
for jj=1:size(Ca_ins,2)
    plot(Vm,I_NaCa(Vm,Ca_ins(jj)),'LineWidth',4); hold on
end
xlabel('membrane potential, mV'); ylabel('current, pA'); set(gca,'FontSize',20);
lgd  = legend(['[Ca]_i = ' num2str(Ca_ins(1)) ' mM'],['[Ca]_i = ' num2str(Ca_ins(2)) ' mM'],...
['[Ca]_i = ' num2str(Ca_ins(3)) ' mM'], ['[Ca]_i = ' num2str(Ca_ins(4)) ' mM'],...
['[Ca]_i = ' num2str(Ca_ins(5)) ' mM'], ['[Ca]_i = ' num2str(Ca_ins(6)) ' mM']);%, ...
% ['[Ca]_i = ' num2str(Ca_ins(7)) ' mM'],['[Ca]_i = ' num2str(Ca_ins(8)) ' mM']);
lgd.Location = 'southeast';
line([-100,50],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-50,200],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
saveas(fig,'figs/gall_susa_1999_model_figs/gall_susa_1999_model_IV_curve_vs_Cain.png');

%I-V vs Na_in
fig=figure(2); clf
V_NaCa = @(Na_in) (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out/Ca_in))*1e3; %mV
I_NaCa = @(Vm,Na_in) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa(Na_in));
Na_ins = [16; 14; 12; 10; 8; 6;];
for jj=1:size(Na_ins,1)
    plot(Vm,I_NaCa(Vm,Na_ins(jj)),'LineWidth',4); hold on
end
xlabel('membrane potential, mV'); ylabel('current, pA'); set(gca,'FontSize',20);
lgd=legend(['[Na]_i = ' num2str(Na_ins(1)) ' mM'],['[Na]_i = ' num2str(Na_ins(2)) ' mM'],...
['[Na]_i = ' num2str(Na_ins(3)) ' mM'], ['[Na]_i = ' num2str(Na_ins(4)) ' mM'],...
['[Na]_i = ' num2str(Na_ins(5)) ' mM'], ['[Na]_i = ' num2str(Na_ins(6)) ' mM']);
lgd.Location = 'southeast';
line([-100,50],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-700,200],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
saveas(fig,'figs/gall_susa_1999_model_figs/gall_susa_1999_model_IV_curve_vs_Nain.png');

%saturation curves
fig=figure(3); clf
V_NaCa = @(Ca_in) (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out./Ca_in))*1e3; %mV
I_NaCa = @(Vm,Ca_in)x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa(Ca_in));
Ca_ins = logspace(-5, -2, 100);
semilogx(Ca_ins,I_NaCa(50,Ca_ins)./max(abs(I_NaCa(50,Ca_ins))),'-','LineWidth',4)
xlabel('[Ca^{2+}]_i, mM'); ylabel('normalized current');
set(gca,'FontSize',20); title('Gall & Susa 1999 Model, [Ca^{2+}]_i saturation curves at 50 mV');
saveas(fig,'figs/gall_susa_1999_model_figs/gall_susa_1999_model_cain_saturation_curve.png');


fig=figure(4); clf
subplot(2,2,3);
V_NaCa = @(Ca_in) (R*T/F)*(3*log(Na_out./Na_in)-log(Ca_out./Ca_in))*1e3; %mV
I_NaCa = @(Vm,Ca_in)x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa(Ca_in));
Ca_ins = logspace(-5, -2, 100);
semilogx(Ca_ins,I_NaCa(50,Ca_ins),'-','LineWidth',4)
xlabel('[Ca^{2+}]_i, mM'); ylabel('current, pA');
set(gca,'FontSize',20);

subplot(2,2,1);
V_NaCa = @(Ca_out) (R*T/F)*(3*log(Na_out./Na_in)-log(Ca_out./Ca_in))*1e3; %mV
I_NaCa = @(Vm,Ca_out)x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa(Ca_out));
Ca_outs = logspace(-3, 2, 100);
semilogx(Ca_outs,I_NaCa(50,Ca_outs),'-','LineWidth',4)
xlabel('[Ca^{2+}]_e, mM'); ylabel('current, pA');
set(gca,'FontSize',20);

subplot(2,2,4);
V_NaCa = @(Na_in) (R*T/F)*(3*log(Na_out./Na_in)-log(Ca_out./Ca_in))*1e3; %mV
I_NaCa = @(Vm,Na_in)x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa(Na_in));
Na_ins = logspace(-1, 4, 100);
semilogx(Na_ins,I_NaCa(50,Na_ins),'-','LineWidth',4)
xlabel('[Na^{+}]_e, mM'); ylabel('current, pA');
set(gca,'FontSize',20);

subplot(2,2,2);
V_NaCa = @(Na_out) (R*T/F)*(3*log(Na_out./Na_in)-log(Ca_out./Ca_in))*1e3; %mV
I_NaCa = @(Vm,Na_out)x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa(Na_out));
Na_outs = logspace(-1, 2, 100);
semilogx(Na_outs,I_NaCa(50,Na_outs),'-','LineWidth',4)
xlabel('[Na^{+}]_e, mM'); ylabel('current, pA');
set(gca,'FontSize',20);

saveas(fig,'figs/gall_susa_1999_model_figs/gall_susa_1999_model_saturation_curves_50mv.png');

fig=figure(5); clf
subplot(2,2,3);
V_NaCa = @(Ca_in) (R*T/F)*(3*log(Na_out./Na_in)-log(Ca_out./Ca_in))*1e3; %mV
I_NaCa = @(Vm,Ca_in)x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa(Ca_in));
Ca_ins = logspace(-5, -2, 100);
semilogx(Ca_ins,I_NaCa(-100,Ca_ins),'-','LineWidth',4)
xlabel('[Ca^{2+}]_i, mM'); ylabel('current, pA');
set(gca,'FontSize',20);

subplot(2,2,1);
V_NaCa = @(Ca_out) (R*T/F)*(3*log(Na_out./Na_in)-log(Ca_out./Ca_in))*1e3; %mV
I_NaCa = @(Vm,Ca_out)x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa(Ca_out));
Ca_outs = logspace(-3, 2, 100);
semilogx(Ca_outs,I_NaCa(-100,Ca_outs),'-','LineWidth',4)
xlabel('[Ca^{2+}]_e, mM'); ylabel('current, pA');
set(gca,'FontSize',20);

subplot(2,2,4);
V_NaCa = @(Na_in) (R*T/F)*(3*log(Na_out./Na_in)-log(Ca_out./Ca_in))*1e3; %mV
I_NaCa = @(Vm,Na_in)x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa(Na_in));
Na_ins = logspace(-1, 4, 100);
semilogx(Na_ins,I_NaCa(-100,Na_ins),'-','LineWidth',4)
xlabel('[Na^{+}]_e, mM'); ylabel('current, pA');
set(gca,'FontSize',20);

subplot(2,2,2);
V_NaCa = @(Na_out) (R*T/F)*(3*log(Na_out./Na_in)-log(Ca_out./Ca_in))*1e3; %mV
I_NaCa = @(Vm,Na_out)x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa(Na_out));
Na_outs = logspace(-1, 2, 100);
semilogx(Na_outs,I_NaCa(-100,Na_outs),'-','LineWidth',4)
xlabel('[Na^{+}]_e, mM'); ylabel('current, pA');
set(gca,'FontSize',20);

saveas(fig,'figs/gall_susa_1999_model_figs/gall_susa_1999_model_saturation_curves_min100mv.png');

