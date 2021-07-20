% EAAT2 - Glumate Transporter phenomenological model
% based on Wadiche et al. 1995b

clear; close all;

%constants
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

%-----------Wadiche et al. 1995b Data
%fit to Wadiche et al., 1995b Fig2B EAAT2 I-V data
[Vm_data, I_data] = csvimport('data/Wadiche_1995b_fig2b_10microM.csv', 'columns', [1, 2] ,'noHeader', true);
%fix positive currents to zero
I_data(I_data>0) = 0;
%ionic concentrations
%intracellular concentrations
%recorded directly from oocytes, so whatever internal concentrations are
%natural - from Gillespie 1983
Na_in = 91; %mM
K_in = 4;  %mM
H_in = 10^(-7); %mM - not sure about this
Glu_in = 0; %mM - not stated
Cl_in = 59; %mM

%extracellular concentrations
%from Wadiche 1995b:
%Recording solution (ND96) contained 96 mM NaCl, 2 mM KCl, 1 mM MgCl2, 1.8 mM CaCl2, and 5 mM HEPES (pH 7.4).
Na_out = 96; %mM
K_out = 2;  %mM
H_out = 10^(-7.4); %mM
Cl_out = 96 + 2 + 1*2 + 1.8*2; %mM
Glu_out = 0.01; %mM

%Vm for I-V curve
Vm = (-140:5:100); %mV

%functions - note RT/F has units volts
%reversal potential for [Cl-]
Vcl = (R*T/(-1*F))*log(Cl_out/Cl_in)*1e3; %mV

%model params
kglu = 1.58e-2; %mM, from wadiche_1995b_model_hillfit
kNa = 46; %mM, from Zerangue & Kavanaugh 1996 Fig 3
%nNa = 2.3, also from Zerangue & Kavanaugh 1996 Fig 3, makes sense since 
%3 Na are transported each cycle of the transporter
kH = 26e-6; %mM, from Zerangue & Kavanaugh 1996 Fig 3
kK = 17; %mM, from Zerangue & Kavanaugh 1996 Fig 3
% gcl = ?; %chloride conductance
%gcl is unknown x(1) below
% tau = 5; %glutamate transporter turnover rate at 0mV, 5 1/s from Wadiche 1995a Fig 7
% N = 4.1e10/(6.022e23); %number of moles of transporter expressed
% tau*N*F is unknown x(2) below
% mu = ?; %Boltzmann factor determining the voltage depen- dence of the transporter turnover rate
% mu is unknown x(3) below

%model
I_EAAT2 = @(x,Vm) (Na_out^2.3./(kNa^2.3+Na_out^2.3))*(H_out./(kH+H_out))*...
    (K_in./(kK+K_in))*(Glu_out./(kglu+Glu_out))*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));


%my fit
my_fit = lsqcurvefit(I_EAAT2,[0, 16,0.02],Vm_data,I_data);

figure(1);clf;
plot(Vm, I_EAAT2(my_fit,Vm), 'b-','LineWidth',3);  hold on;
xlabel('membrane potential, mV'); ylabel('current, pA');
set(gca,'FontSize',20); title('EAAT2 Data (Wadiche et al. 1995b Fig 2B)');
xlim([-125,140]); ylim([-100,10]);
plot(Vm_data, I_data, 'bo','MarkerSize',10,'LineWidth',2,'MarkerFaceColor','b')
legend('model', 'Wadiche et al 1995b Fig1D data')
%axes lines
line([-165,140],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');


%--- testing fit vs. more [Glu]out
%30 microM
[Vm_data2, I_data2] = csvimport('data/Wadiche_1995b_fig2b_30microM.csv', 'columns', [1, 2] ,'noHeader', true);
%fix positive currents to zero
I_data2(I_data2>0) = 0;
Glu_out = 0.03; %mM
%remake model with new [Glu]out
I_EAAT2 = @(x,Vm) (Na_out^2.3./(kNa^2.3+Na_out^2.3))*(H_out./(kH+H_out))*...
    (K_in./(kK+K_in))*(Glu_out./(kglu+Glu_out))*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));

plot(Vm, I_EAAT2(my_fit,Vm), 'r-','LineWidth',3); 
plot(Vm_data2, I_data2, 'r^','MarkerSize',10,'LineWidth',2)

%100 microM
[Vm_data3, I_data3] = csvimport('data/Wadiche_1995b_fig2b_100microM.csv', 'columns', [1, 2] ,'noHeader', true);
%fix positive currents to zero
I_data3(I_data3>0) = 0;
Glu_out = 0.1; %mM
%remake model with new [Glu]out
I_EAAT2 = @(x,Vm) (Na_out^2.3./(kNa^2.3+Na_out^2.3))*(H_out./(kH+H_out))*...
    (K_in./(kK+K_in))*(Glu_out./(kglu+Glu_out))*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));

plot(Vm, I_EAAT2(my_fit,Vm), 'm-','LineWidth',3);  
plot(Vm_data3, I_data3, 'md','MarkerSize',10,'LineWidth',2)

%300 microM
[Vm_data4, I_data4] = csvimport('data/Wadiche_1995b_fig2b_300microM.csv', 'columns', [1, 2] ,'noHeader', true);
%fix positive currents to zero
I_data4(I_data4>0) = 0;
Glu_out = 0.3; %mM
%remake model with new [Glu]out
I_EAAT2 = @(x,Vm) (Na_out^2.3./(kNa^2.3+Na_out^2.3))*(H_out./(kH+H_out))*...
    (K_in./(kK+K_in))*(Glu_out./(kglu+Glu_out))*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));

plot(Vm, I_EAAT2(my_fit,Vm), 'c-','LineWidth',3);  
plot(Vm_data4, I_data4, 'cs','MarkerSize',10,'LineWidth',2)

%1 mM
[Vm_data5, I_data5] = csvimport('data/Wadiche_1995b_fig2b_1mM.csv', 'columns', [1, 2] ,'noHeader', true);
%fix positive currents to zero
I_data5(I_data5>0) = 0;
Glu_out = 1; %mM
%remake model with new [Glu]out
I_EAAT2 = @(x,Vm) (Na_out^2.3./(kNa^2.3+Na_out^2.3))*(H_out./(kH+H_out))*...
    (K_in./(kK+K_in))*(Glu_out./(kglu+Glu_out))*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));

plot(Vm, I_EAAT2(my_fit,Vm), 'k-','LineWidth',3);  
plot(Vm_data5, I_data5, 'ko','MarkerSize',10,'LineWidth',2)
legend('model', '0.01 mM Glu', 'model', '0.03 mM Glu',...
    'model', '0.1 mM Glu', 'model', '0.3 mM Glu', 'model', '1 mM Glu');


%--- testing vs. different ion concentration
Vm = -80; %mV
Na_out = 96; %mM
H_out = 10^(-7.4); %mM
K_in = 4; %mM
Glu_out = linspace(0,.300);
I_EAAT2 = @(x,Vm, Glu_out) (Na_out.^2.3./(kNa^2.3+Na_out.^2.3)).*(H_out./(kH+H_out)).*...
    (K_in./(kK+K_in)).*(Glu_out./(kglu+Glu_out)).*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));
figure(3); clf;
subplot(2,2,1);
plot(Glu_out*1e3, ...
    abs(I_EAAT2(my_fit,Vm,Glu_out))./max(abs(I_EAAT2(my_fit,Vm,Glu_out))),...
    'k-','LineWidth',3);  
xlabel('[Glu]_o \muM'); ylabel('normalized current'); set(gca,'FontSize',20);

figure(4); clf;
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,Glu_out(2)),'-','LineWidth',3); hold on;  
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,Glu_out(4)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,Glu_out(10)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,Glu_out(20)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,Glu_out(40)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,Glu_out(100)),'-','LineWidth',3);
xlim([-125,140]); ylim([-100,10]);
xlabel('membrane potential, mV'); ylabel('current, pA'); set(gca,'FontSize',20);
legend(['[Glu]_o = ' num2str(Glu_out(2)) ' mM'],['[Glu]_o = ' num2str(Glu_out(20)) ' mM'],...
    ['[Glu]_o = ' num2str(Glu_out(40)) ' mM'], ['[Glu]_o = ' num2str(Glu_out(60)) ' mM'],...
    ['[Glu]_o = ' num2str(Glu_out(80)) ' mM'], ['[Glu]_o = ' num2str(Glu_out(100)) ' mM']);
line([-165,140],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');


%-- Na_out
Vm = -80; %mV
Na_out = linspace(0,150);
H_out = 10^(-7.4); %mM
K_in = 4; %mM
Glu_out = 0.01; %mM
I_EAAT2 = @(x,Vm, Na_out) (Na_out.^2.3./(kNa^2.3+Na_out.^2.3)).*(H_out./(kH+H_out)).*...
    (K_in./(kK+K_in)).*(Glu_out./(kglu+Glu_out))*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));
figure(3);
subplot(2,2,2);
plot(Na_out, ...
    abs(I_EAAT2(my_fit,Vm,Na_out))./max(abs(I_EAAT2(my_fit,Vm,Na_out))),...
    'k-','LineWidth',3);  
xlabel('[Na]_o mM'); ylabel('normalized current'); set(gca,'FontSize',20);
figure(6); clf;
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,Na_out(2)),'-','LineWidth',3); hold on;  
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,Na_out(20)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,Na_out(40)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,Na_out(60)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,Na_out(80)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,Na_out(100)),'-','LineWidth',3);
xlim([-125,140]); ylim([-100,10]);
xlabel('membrane potential, mV'); ylabel('current, pA'); set(gca,'FontSize',20);
legend(['[Na]_o = ' num2str(Na_out(2)) ' mM'],['[Na]_o = ' num2str(Na_out(20)) ' mM'],...
    ['[Na]_o = ' num2str(Na_out(40)) ' mM'], ['[Na]_o = ' num2str(Na_out(60)) ' mM'],...
    ['[Na]_o = ' num2str(Na_out(80)) ' mM'], ['[Na]_o = ' num2str(Na_out(100)) ' mM']);
line([-165,140],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

%-- H_out
Vm = -80; %mV
Na_out = 96; %mM
H_out = logspace(-7.4,-4,100); %mM
K_in = 4; %mM
Glu_out = 0.01; %mM
I_EAAT2 = @(x,Vm, H_out) (Na_out.^2.3./(kNa^2.3+Na_out^2.3)).*(H_out./(kH+H_out)).*...
    (K_in./(kK+K_in)).*(Glu_out./(kglu+Glu_out)).*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));
figure(3); 
subplot(2,2,3);
plot(H_out*1e6, ...
    abs(I_EAAT2(my_fit,Vm,H_out))./max(abs(I_EAAT2(my_fit,Vm,H_out))),...
    'k-','LineWidth',3);  
xlabel('[H]_o nM'); ylabel('normalized current'); set(gca,'FontSize',20);
figure(8); clf;
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,H_out(2)),'-','LineWidth',3); hold on;  
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,H_out(4)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,H_out(10)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,H_out(20)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,H_out(40)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,H_out(60)),'-','LineWidth',3);
xlim([-125,140]); ylim([-100,10]);
xlabel('membrane potential, mV'); ylabel('current, pA'); set(gca,'FontSize',20);
legend(['[H]_o = ' num2str(H_out(2)) ' mM'],['[H]_o = ' num2str(H_out(20)) ' mM'],...
    ['[H]_o = ' num2str(H_out(40)) ' mM'], ['[H]_o = ' num2str(H_out(60)) ' mM'],...
    ['[H]_o = ' num2str(H_out(80)) ' mM'],['[H]_o = ' num2str(H_out(100)) ' mM']);
line([-165,140],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

%-- K_in
Vm = -80; %mV
Na_out = 96; %mM
H_out = 10^(-7.4); %mM
K_in = linspace(0,70); %mM
Glu_out = 0.01; %mM
I_EAAT2 = @(x,Vm, K_in) (Na_out^2.3./(kNa^2.3+Na_out^2.3)).*(H_out./(kH+H_out)).*...
    (K_in./(kK+K_in)).*(Glu_out./(kglu+Glu_out)).*((Vm-Vcl)*x(1) - x(2)*exp(-x(3)*Vm));
figure(3);
subplot(2,2,4);
plot(K_in, ...
    abs(I_EAAT2(my_fit,Vm,K_in))./max(abs(I_EAAT2(my_fit,Vm,K_in))),...
    'k-','LineWidth',3);  
xlabel('[K]_i mM'); ylabel('normalized current'); set(gca,'FontSize',20);
figure(10); clf;
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,K_in(2)),'-','LineWidth',3); hold on;  
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,K_in(10)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,K_in(20)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,K_in(40)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,K_in(60)),'-','LineWidth',3);
plot(-140:5:100, I_EAAT2(my_fit,-140:5:100,K_in(100)),'-','LineWidth',3);
xlim([-125,140]); ylim([-100,10]);
xlabel('membrane potential, mV'); ylabel('current, pA'); set(gca,'FontSize',20);
legend(['[K]_i = ' num2str(K_in(2)) ' mM'],['[K]_i = ' num2str(K_in(20)) ' mM'],...
    ['[K]_i = ' num2str(K_in(40)) ' mM'], ['[K]_i = ' num2str(K_in(60)) ' mM'],...
    ['[K]_i = ' num2str(K_in(80)) ' mM'], ['[K]_i = ' num2str(K_in(100)) ' mM']);
line([-165,140],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k','HandleVisibility','off');