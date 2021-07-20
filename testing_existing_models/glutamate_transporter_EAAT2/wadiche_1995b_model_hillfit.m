% EAAT2 - Glumate Transporter model
% based on Wadiche et al. 1995b
% Hill fit for current

clear; close all;

%-----------Wadiche et al. 1995b Data
%fit to Wadiche et al., 1995b Fig2B EAAT2 I-V data
%10 microM
[Vm_data, I_data] = csvimport('data/Wadiche_1995b_fig2b_10microM.csv', 'columns', [1, 2] ,'noHeader', true);
%30 microM
[Vm_data2, I_data2] = csvimport('data/Wadiche_1995b_fig2b_30microM.csv', 'columns', [1, 2] ,'noHeader', true);
%100 microM
[Vm_data3, I_data3] = csvimport('data/Wadiche_1995b_fig2b_100microM.csv', 'columns', [1, 2] ,'noHeader', true);
%300 microM
[Vm_data4, I_data4] = csvimport('data/Wadiche_1995b_fig2b_300microM.csv', 'columns', [1, 2] ,'noHeader', true);
%1 mM
[Vm_data5, I_data5] = csvimport('data/Wadiche_1995b_fig2b_1mM.csv', 'columns', [1, 2] ,'noHeader', true);

Hill_fit_Is = [I_data(1); I_data2(1); I_data3(1); I_data4(1); I_data5(1);];
Hill_fit_cs = [0.01; 0.03; 0.1; 0.3; 1;]; %mM

I_norm = @(x,Glu) Glu./(x+Glu)
hill_fit = lsqcurvefit(I_norm,1.5e-2,Hill_fit_cs,Hill_fit_Is./min(Hill_fit_Is));
figure(2); 
plot(Hill_fit_cs, Hill_fit_Is./min(Hill_fit_Is)); hold on;
plot(Hill_fit_cs, I_norm(hill_fit,Hill_fit_cs))
hill_fit