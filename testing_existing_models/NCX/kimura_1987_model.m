% NCX - Sodium-Calcium exchanger model
% from Kimura et al. 1987, based on DiFrancesco & Noble, 1985
% still need to fit to astrocyte data


%constants
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 308.15; %K, absolute temperature (35 C)

%cellular concentrations
Na_in = 20; %mM
Ca_in = 67*1e-6; %mM
Na_out = 140; %mM
Ca_out = 0.1; %mM

%functions - note RT/F has units volts
%reversal potential for NCX
V_NaCa = (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out/Ca_in))*1e3; %mV
V_NaCa


%hill functions parameters
x = [2.07e-5; 0.38;]; %muA muF /mM^4, partition coefficient (symmetric reversal potential)
I_NaCa = @(Vm) x(1)*((Na_in.^3).*(Ca_out).*(exp(x(2)*Vm*1e-3*F/(R*T))) ...
   - (Na_out.^3).*(Ca_in).*(exp(-(1-x(2))*Vm*1e-3*F/(R*T)))) ;

%Vm for I-V curve
Vm = (-100:5:50); %mV

fig = figure(1); clf
I_NaCa = @(Vm, Ca_in) x(1)*((Na_in.^3).*(Ca_out).*(exp(x(2)*Vm*1e-3*F/(R*T))) ...
   - (Na_out.^3).*(Ca_in).*(exp(-(1-x(2))*Vm*1e-3*F/(R*T)))) ;
Ca_ins = [600; 500; 400; 300; 200; 100;].*1e-6;
for jj=1:size(Ca_ins,1)
    plot(Vm,I_NaCa(Vm,Ca_ins(jj)),'LineWidth',4); hold on
end
% plot(Vm,I_NaCa(Vm),'LineWidth',4);  hold on
% plot([V_NaCa],[0],'o','LineWidth',4);
% end
% xlim([-125,140]); ylim([-200,200]);
xlabel('membrane potential, mV'); ylabel('current, \muA'); set(gca,'FontSize',20);
lgd= legend(['[Ca]_i = ' num2str(Ca_ins(1)) ' mM'],['[Ca]_i = ' num2str(Ca_ins(2)) ' mM'],...
['[Ca]_i = ' num2str(Ca_ins(3)) ' mM'], ['[Ca]_i = ' num2str(Ca_ins(4)) ' mM'],...
['[Ca]_i = ' num2str(Ca_ins(5)) ' mM'], ['[Ca]_i = ' num2str(Ca_ins(6)) ' mM']);
lgd.Location= 'southeast';
line([-100,50],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-.35,0.05],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

saveas(fig,'figs/kimura_etal_1987_model_figs/kimura_etal_1987_model_IV_curve_vsCain.png');

fig = figure(2); clf
I_NaCa = @(Vm, Na_in) x(1)*((Na_in.^3).*(Ca_out).*(exp(x(2)*Vm*1e-3*F/(R*T))) ...
   - (Na_out.^3).*(Ca_in).*(exp(-(1-x(2))*Vm*1e-3*F/(R*T)))) ;
Na_ins = [16; 14; 12; 10; 8; 6;];
for jj=1:size(Na_ins,1)
    plot(Vm,I_NaCa(Vm,Na_ins(jj)),'LineWidth',4); hold on
end
xlabel('membrane potential, mV'); ylabel('current, \muA'); set(gca,'FontSize',20);
lgd= legend(['[Na]_i = ' num2str(Na_ins(1)) ' mM'],['[Na]_i = ' num2str(Na_ins(2)) ' mM'],...
['[Na]_i = ' num2str(Na_ins(3)) ' mM'], ['[Na]_i = ' num2str(Na_ins(4)) ' mM'],...
['[Na]_i = ' num2str(Na_ins(5)) ' mM'], ['[Na]_i = ' num2str(Na_ins(6)) ' mM']);
lgd.Location= 'southeast';
line([-100,50],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-.05,0.05],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

saveas(fig,'figs/kimura_etal_1987_model_figs/kimura_etal_1987_model_IV_curve_vsNain.png');


%saturation curves
fig = figure(3);clf
subplot(2,2,1);
I_NaCa = @(Vm, Ca_out) x(1)*((Na_in.^3).*(Ca_out).*(exp(x(2)*Vm*1e-3*F/(R*T))) ...
   - (Na_out.^3).*(Ca_in).*(exp(-(1-x(2))*Vm*1e-3*F/(R*T)))) ;

Ca_outs = logspace(-1,2.15,100);
semilogx(Ca_outs,abs(I_NaCa(-100,Ca_outs)),'LineWidth',4); hold on
% xlim([-125,140]); ylim([-200,200]);
xlabel('[Ca^{2+}]_e, mM'); ylabel('current, \muA'); set(gca,'FontSize',20);

subplot(2,2,3);
I_NaCa = @(Vm, Ca_in) x(1)*((Na_in.^3).*(Ca_out).*(exp(x(2)*Vm*1e-3*F/(R*T))) ...
   - (Na_out.^3).*(Ca_in).*(exp(-(1-x(2))*Vm*1e-3*F/(R*T)))) ;

Ca_ins = logspace(-1,2.15,100);
semilogx(Ca_outs,abs(I_NaCa(-100,Ca_ins)),'LineWidth',4); hold on
% xlim([-125,140]); ylim([-200,200]);
xlabel('[Ca^{2+}]_i, mM'); ylabel('current, \muA'); set(gca,'FontSize',20);

subplot(2,2,4);
I_NaCa = @(Vm, Na_in) x(1)*((Na_in.^3).*(Ca_out).*(exp(x(2)*Vm*1e-3*F/(R*T))) ...
   - (Na_out.^3).*(Ca_in).*(exp(-(1-x(2))*Vm*1e-3*F/(R*T)))) ;

Na_ins = logspace(-1,2.15,100);
semilogx(Ca_outs,abs(I_NaCa(-100,Na_ins)),'LineWidth',4); hold on
% xlim([-125,140]); ylim([-200,200]);
xlabel('[Na^{+}]_i, mM'); ylabel('current, \muA'); set(gca,'FontSize',20);



%matching fig 8 - Na_out saturation curve
subplot(2,2,2);
%cellular concentrations
Na_in = 0; %mM
Ca_in = 430*1e-6; %mM
Na_out = 140; %mM
Ca_out = 1; %mM
I_NaCa = @(Vm, Na_out) x(1)*((Na_in.^3).*(Ca_out).*(exp(x(2)*Vm*1e-3*F/(R*T))) ...
   - (Na_out.^3).*(Ca_in).*(exp(-(1-x(2))*Vm*1e-3*F/(R*T)))) ;

Na_outs = logspace(1.24,2.15,100);

semilogx(Na_outs,abs(I_NaCa(-100,Na_outs)),'LineWidth',4); hold on
% xlim([-125,140]); ylim([-200,200]);
xlabel('[Na^+]_e, mM'); ylabel('current, \muA'); set(gca,'FontSize',20);

saveas(fig,'figs/kimura_etal_1987_model_figs/kimura_etal_1987_model_saturationcurves.png');



