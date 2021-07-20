% EAAT2 - Glumate Transporter model
% based on I-V data

clear; close all; clc;
%Kir
%constants
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

%-----------Levy et al., 1998 Data
%ionic concentrations - data from Levy et al. 1998
%intracellular concentrations
Na_in = 12; %mM
K_in = 140;  %mM
H_in = 10^(-7); %mM
Glu_in = 0.2; %mM
%extracellular concentrations
Na_out = 142; %mM
K_out = 2.5;  %mM
H_out = 10^(-7.4); %mM
Glu_out = 0.1; %mM

%flanagan et al. 2018 model params
alpha = 1.9767*10^(-5); %A m^-2
beta = .0292; %mV^-1

%Vm for I-V curve
Vm = (-125:5:140); %mV

%functions - note RT/F has units volts
V_rev = (R*T/(2*F))*log((Na_out/Na_in)^3*(K_in/K_out)*(H_out/H_in)*(Glu_out/Glu_in))*1e3; %mV

%model from Flanagan et al. 2018
I_EAAT2_flan = @(x,Vm) -x(1).*exp(x(2).*(Vm-V_rev));
%my model
kglu = 1.59e-2;
nglu = 1.1141;
% kglu = 0.0730; %halfway between each fit
% nglu = 1;

I_EAAT2 = @(x,Vm) x(1).*exp(x(2).*(Vm-x(3)-V_rev)).*(Glu_out.^nglu./(kglu.^nglu+Glu_out.^nglu));

%fit to Levy et al., 1998 Fig1 (B) Glutate-transporter GLT-1 (EAAT2 analog) I-V data
[Vm_data, I_data] = csvimport('data/levy_etal_1998_fig1b.csv', 'columns', [1, 2] ,'noHeader', true);
%my fit
my_fit = lsqcurvefit(I_EAAT2,[0.0511, -0.0134,0],Vm_data,I_data);
%add (V_rev, 0) to dataset
Vm_data(end+1) = V_rev;
I_data(end+1) = 0;
%flagan fit
x0 = [alpha,beta];
flanagan_fit = lsqcurvefit(I_EAAT2_flan,x0,Vm_data,I_data);

figure(1);clf;
plot(Vm, I_EAAT2_flan(flanagan_fit,Vm), 'b-','LineWidth',3); hold on
plot(Vm, I_EAAT2(my_fit,Vm), 'r-','LineWidth',3); 
xlabel('membrane potential, mV'); ylabel('current, pA');
set(gca,'FontSize',20);
xlim([-125,140]); ylim([-100,0]);
plot(Vm_data, I_data, 'm^','MarkerSize',10,'LineWidth',2)
legend('Flanagan et al., 2018 model', 'my model', 'Levy et al., 1998 data')
%axes lines
line([-165,140],[0,0],'LineStyle','--', 'Color', 'k');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k');


%--------------- Testing models vs. Extracellular Glutamate concentration ----
[Glu_out_data, I_Glu_data] = csvimport('data/levy_etal_1998_fig2b.csv', 'columns', [1, 2] ,'noHeader', true);
Glu_out_data = Glu_out_data*1e-3; %convert from microM to mM
Glu_out = logspace(log10(min(Glu_out_data)),log10(max(Glu_out_data)),20);
Vm = -110;
%new V_reversal, glutamate-dependent
V_rev = @(Glu) (R*T/(2*F)).*log((Na_out/Na_in)^3*(K_in/K_out)*(H_out/H_in).*(Glu./Glu_in))*1e3; %mV
%remake Flanagan functions with new V_rev since it depends on it
I_EAAT2_flan = @(x,Vm,Glu) -x(1).*exp(x(2).*(Vm-V_rev(Glu)));
%remake my model function with new Glu_out since it depends on it
I_EAAT2 = @(x,Vm,Glu) x(1).*exp(x(2).*(Vm-x(3)-V_rev(Glu)))...
                .*(Glu.^nglu./(kglu.^nglu+Glu.^nglu));

%plot
figure(3);clf;
semilogx(Glu_out, -I_EAAT2_flan(flanagan_fit,Vm,Glu_out)./max(abs(I_EAAT2_flan(flanagan_fit,Vm,Glu_out))),...
    'b-','LineWidth',3); hold on
semilogx(Glu_out, -I_EAAT2(my_fit,Vm,Glu_out)./max(abs(I_EAAT2(my_fit,Vm,Glu_out))), 'r-','LineWidth',3);
xlabel('extracellular Glutamate concentration, mM'); ylabel('current, normalized');
set(gca,'FontSize',20);
semilogx(Glu_out_data, I_Glu_data, 'm^','MarkerSize',10,'LineWidth',2)
legend('Flanagan et al., 2018 model', 'my model', 'Levy et al., 1998 data')

return
%------------Bergles and Jahr, 1997 Data
%fit to Bergles and Jahr, 1997 Fig7(B) Glutate-transporter EAAT2 I-V data
[Vm_data2, I_data2] = csvimport('data/bergles_jahr_1997_fig7b.csv', 'columns', [1, 2] ,'noHeader', true);
%scale current data to other paper's data
I_data2 = I_data2.*max(abs(I_data));

%Vm for I-V curve
Vm = (-125:5:140); %mV

%intracellular concentrations
Na_in = 0.01; %mM %total guess based on other data, not reported in paper
K_in = 130;  %mM
H_in = 10^(-7.2); %mM
Glu_in = 1; %mM %total guess based on other data, not reported in paper
%extracellular concentrations
Na_out = 135; %mm
K_out = 5.4;  %mM %typo in Flanagan's SI, they say they used 4.5mM but Bergles & Jahr report 5.4mM
H_out = 10^(-7.2); %mM
Glu_out = 10; %mM %also a guess based on the paper mentioning 10mM glutamate used to activate the channel

%back to single V_reversal
V_rev = (R*T/(2*F))*log((Na_out/Na_in)^3*(K_in/K_out)*(H_out/H_in)*(Glu_out/Glu_in))*1e3; %mV
%model from Flanagan et al. 2018
I_EAAT2_flan = @(x,Vm) -x(1).*exp(x(2).*(Vm-V_rev));
%my model
% kglu = 1.3e-1;
% nglu = 0.6479;

I_EAAT2 = @(x,Vm) x(1).*(Vm-x(3)).*exp(x(2).*(Vm-x(3))).*(Glu_out.^nglu./(kglu.^nglu+Glu_out.^nglu));
%my fit
my_fit2 = lsqcurvefit(I_EAAT2,my_fit,Vm_data2,I_data2);

%add (V_rev, 0) to dataset
Vm_data2(end+1) = V_rev;
I_data2(end+1) = 0;
%flanagan fit
flanagan_fit2 = lsqcurvefit(I_EAAT2_flan,flanagan_fit,Vm_data2,I_data2);

figure(2);clf;
plot(Vm, I_EAAT2_flan(flanagan_fit2,Vm), 'b-','LineWidth',3); hold on
plot(Vm, I_EAAT2(my_fit2,Vm), 'r-','LineWidth',3); 
xlabel('membrane potential, mV'); ylabel('current, pA');
set(gca,'FontSize',20);
xlim([-125,140]); ylim([-100,0]);
plot(Vm_data2, I_data2, 'm^','MarkerSize',10,'LineWidth',2)
legend('Flanagan et al., 2018 model', 'my model', 'Bergles and Jahr, 1997 data')
%axes lines
line([-165,140],[0,0],'LineStyle','--', 'Color', 'k');
line([0,0],[-150,100],'LineStyle','--', 'Color', 'k');

%---- Testing models vs. Extracellular Glutamate w/ Bergles & Jahr data
[Glu_out_data, I_Glu_data] = csvimport('data/bergles_jahr_1997_fig8c.csv', 'columns', [1, 2] ,'noHeader', true);
[Glu_out_data2, I_Glu_data2] = csvimport('data/bergles_jahr_1997_fig8c_peak.csv', 'columns', [1, 2] ,'noHeader', true);
Glu_out = logspace(log10(min(Glu_out_data)),log10(max(Glu_out_data)),20);
Vm = -100;

%new V_reversal, glutamate-dependent
V_rev = (R*T/(2*F)).*log((Na_out/Na_in)^3*(K_in/K_out)*(H_out/H_in).*(Glu_out./Glu_in))*1e3; %mV
%remake Flanagan functions with new V_rev since it depends on it
I_EAAT2_flan = @(x,Vm,V_r) -x(1).*exp(x(2).*(Vm-V_r));
%remake my model function with new Glu_out since it depends on it
I_EAAT2 = @(x,Vm,Glu) x(1).*(Vm-x(3)).*exp(x(2).*(Vm-x(3))).*(Glu.^nglu./(kglu.^nglu+Glu.^nglu));
%plot
figure(4);clf;
semilogx(Glu_out, -I_EAAT2_flan(flanagan_fit2,Vm,V_rev)./max(abs(I_EAAT2_flan(flanagan_fit2,Vm,V_rev))),...
    'b-','LineWidth',3); hold on
semilogx(Glu_out, -I_EAAT2(my_fit2,Vm,Glu_out)./max(abs(I_EAAT2(my_fit2,Vm,Glu_out))), 'r-','LineWidth',3);
xlabel('extracellular Glutamate concentration, mM'); ylabel('current, normalized');
set(gca,'FontSize',20);
semilogx(Glu_out_data2, I_Glu_data2, 'm^','MarkerSize',10,'LineWidth',2)
legend('Flanagan et al., 2018 model', 'my model', ...
   'Bergles and Jahr, 1997 data, peak')


%---- Testing models vs. Extracellular Glutamate w/ Bergles & Jahr data and
%model fit to Levy et al I-V curve
[Glu_out_data, I_Glu_data] = csvimport('data/bergles_jahr_1997_fig8c.csv', 'columns', [1, 2] ,'noHeader', true);
[Glu_out_data2, I_Glu_data2] = csvimport('data/bergles_jahr_1997_fig8c_peak.csv', 'columns', [1, 2] ,'noHeader', true);
Glu_out = logspace(log10(min(Glu_out_data)),log10(max(Glu_out_data)),20);
Vm = -100;

%new V_reversal, glutamate-dependent
V_rev = (R*T/(2*F)).*log((Na_out/Na_in)^3*(K_in/K_out)*(H_out/H_in).*(Glu_out./Glu_in))*1e3; %mV
%remake Flanagan functions with new V_rev since it depends on it
I_EAAT2_flan = @(x,Vm,V_r) -x(1).*exp(x(2).*(Vm-V_r));
%remake my model function with new Glu_out since it depends on it
I_EAAT2 = @(x,Vm,Glu) x(1).*(Vm-x(3)).*exp(x(2).*(Vm-x(3))).*(Glu.^nglu./(kglu.^nglu+Glu.^nglu));
%plot
figure(5);clf;
semilogx(Glu_out, -I_EAAT2_flan(flanagan_fit,Vm,V_rev)./max(abs(I_EAAT2_flan(flanagan_fit,Vm,V_rev))),...
    'b-','LineWidth',3); hold on
semilogx(Glu_out, -I_EAAT2(my_fit,Vm,Glu_out)./max(abs(I_EAAT2(my_fit,Vm,Glu_out))), 'r-','LineWidth',3);
xlabel('extracellular Glutamate concentration, mM'); ylabel('current, normalized');
set(gca,'FontSize',20);
semilogx(Glu_out_data2, I_Glu_data2, 'm^','MarkerSize',10,'LineWidth',2)
legend('Flanagan et al., 2018 model', 'my model', ...
     'Bergles and Jahr, 1997 data, peak')

%try to do it at same axis/scale as Flanagan SI
%get it to match exactly at 1mM
Glu_out_for_norm_factor = 1;
norm_factor = I_EAAT2(my_fit,Vm,Glu_out_for_norm_factor)./I_Glu_data(3);
V_rev_for_norm_factor = ...
    (R*T/(2*F)).*log((Na_out/Na_in)^3*(K_in/K_out)*(H_out/H_in).*(Glu_out_for_norm_factor./Glu_in))*1e3; %mV
norm_factor2 = I_EAAT2(my_fit,Vm,Glu_out_for_norm_factor)...
    ./I_EAAT2_flan(flanagan_fit,Vm,V_rev_for_norm_factor);

figure(8);clf;
plot(Glu_out, I_EAAT2_flan(flanagan_fit,Vm,V_rev).*norm_factor2,...
    'b-','LineWidth',3); hold on
plot(Glu_out, I_EAAT2(my_fit,Vm,Glu_out), 'r-','LineWidth',3);
xlabel('extracellular Glutamate concentration, mM'); ylabel('current, set to match my model at [Glu]_e=1mM');
set(gca,'FontSize',20);
plot(Glu_out_data2, I_Glu_data2*norm_factor, 'm^','MarkerSize',10,'LineWidth',2)
legend('Flanagan et al., 2018 model', 'my model', ...
    'Bergles and Jahr, 1997 data')


