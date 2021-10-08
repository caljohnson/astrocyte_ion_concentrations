%Matsuoka & Hilgemann 1992 data test
%testing each model against the data and I-V curves

%Fig 9A - I-V curves
%recording at 34-36 deg C (described as exactly the same method as
%Hilgemann et al. 1992)
%recording of an outward exchange current at various membrane potentials
%and various internal sodium concentrations
%Na_in = 100; %mM varied here
Ca_in = 1e-3; %mM (1 microM)
Na_out = 150; %mM
Ca_out = 2; %mM
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 308.15; %K, 35 deg C

%Vm for I-V curve
Vm = (-150:5:100); %mV

%data from Fig9A - first points
[Vm_data, I_data] = csvimport('data/matsuoka_hilgemann_1992_fig9a_100.csv',...
                    'columns', [1, 2] ,'noHeader', true);
Na_in = 100;%mM                
figure(1); clf;
plot(Vm_data,I_data,'o','LineWidth',4,'MarkerSize',10); hold on
xlabel('V_m (mV)'); ylabel('current, pA'); set(gca,'FontSize',20);

%testing egger & niggli 1999 model
%x = [0.2; 0.5;]; %A l^4 /mol^4, partition coefficient (symmetric reversal potential)
I_NaCa_egger = @(x,Vm) x*((Na_in.*1e-3).^3.*(Ca_out.*1e-3).*(exp(0.5*Vm*1e-3*F/(R*T))) ...
   - ((Na_out.*1e-3).^3.*(Ca_in.*1e-3).*(exp(-0.5*Vm*1e-3*F/(R*T)))));
egger_fit = lsqcurvefit(I_NaCa_egger,2e4,Vm_data,I_data);
plot(Vm,I_NaCa_egger(egger_fit,Vm),'-','LineWidth',4);


%testing kimura et al 1987 model
% x = [2.07e-5; 0.38;]; %muA muF /mM^4, partition coefficient (symmetric reversal potential)
I_NaCa_kimura = @(x,Vm) x*((Na_in.^3).*(Ca_out).*(exp(0.38*Vm*1e-3*F/(R*T))) ...
   - (Na_out.^3).*(Ca_in).*(exp(-(1-0.38)*Vm*1e-3*F/(R*T)))) ;
kimura_fit = lsqcurvefit(I_NaCa_kimura,2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_kimura(kimura_fit,Vm),'-.','LineWidth',4);

%testing gall & susa 1999 model
%hill params for gall & susa
% nH = 1;
% K_H = 1.5*1e-3; %1.5 microM
nH = 3.7;
K_H = 0.6*1e-3;%mM
V_NaCa = (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out/Ca_in))*1e3; %mV
I_NaCa_gall_susa = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa);
gall_susa_fit = lsqcurvefit(I_NaCa_gall_susa,-0.07,Vm_data,I_data);
plot(Vm,I_NaCa_gall_susa(gall_susa_fit,Vm),'--','LineWidth',4);

%testing flanagan et al 2018 model
I_NaCa_flanagan = @(x,Vm) x*((Na_in.^3)./(Na_out.^3).*(exp(0.5*Vm*1e-3*F/(R*T))) ...
   - (Ca_in)./(Ca_out).*(exp(0.5*Vm*1e-3*F/(R*T)))) ;
flanagan_fit = lsqcurvefit(I_NaCa_flanagan,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_flanagan(flanagan_fit,Vm),':','LineWidth',4);

%testing my model
nH2 = 1.6;
kH2 = 20;
I_NaCa_my1 = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH) ...
    .*(Na_in.^nH2)./(Na_in.^nH2 + kH2^nH2).*(Vm-V_NaCa);
Cainhill_and_Nainhill = lsqcurvefit(I_NaCa_my1,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_my1(Cainhill_and_Nainhill,Vm),':','LineWidth',4);

I_NaCa_linear = @(x,Vm) x.*(Vm-V_NaCa);
linear_fit = lsqcurvefit(I_NaCa_linear,2,Vm_data,I_data);
plot(Vm,I_NaCa_linear(linear_fit,Vm),':','LineWidth',4);


lgd = legend('Matsuoka & Hilgemann 1992 Fig9A data','Egger & Niggli 1999 model fit',...
    'Kimura et al. 1987 model fit', 'Gall & Susa 1999 model fit',...
    'Flanagan et al. 2018 model fit', 'Ca_i Hill + Na_i Hill', 'Linear');
lgd.Location = 'northwest';

line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-50,50],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

title('[Na^+]_i = 100 mM, [Ca^{2+}]_i = 1 \muM');
xlim([-150,100]); ylim([-50,100])


%data from Fig9A - 50 mM Nain
[V_data, I_data] = csvimport('data/matsuoka_hilgemann_1992_fig9a_50.csv',...
                    'columns', [1, 2] ,'noHeader', true);              
Na_in = 50;%mM          
figure(2); clf;
plot(Vm_data,I_data,'o','LineWidth',4,'MarkerSize',10); hold on
xlabel('V_m (mV)'); ylabel('current, pA'); set(gca,'FontSize',20);

%testing egger & niggli 1999 model
%x = [0.2; 0.5;]; %A l^4 /mol^4, partition coefficient (symmetric reversal potential)
I_NaCa_egger = @(x,Vm) x*((Na_in.*1e-3).^3.*(Ca_out.*1e-3).*(exp(0.5*Vm*1e-3*F/(R*T))) ...
   - ((Na_out.*1e-3).^3.*(Ca_in.*1e-3).*(exp(-0.5*Vm*1e-3*F/(R*T)))));
% egger_fit = lsqcurvefit(I_NaCa_egger,2e3,Vm_data,I_data);
plot(Vm,I_NaCa_egger(egger_fit,Vm),'-','LineWidth',4);


%testing kimura et al 1987 model
% x = [2.07e-5; 0.38;]; %muA muF /mM^4, partition coefficient (symmetric reversal potential)
I_NaCa_kimura = @(x,Vm) x*((Na_in.^3).*(Ca_out).*(exp(0.38*Vm*1e-3*F/(R*T))) ...
   - (Na_out.^3).*(Ca_in).*(exp(-(1-0.38)*Vm*1e-3*F/(R*T)))) ;
% kimura_fit = lsqcurvefit(I_NaCa_kimura,2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_kimura(kimura_fit,Vm),'-.','LineWidth',4);

%testing gall & susa 1999 model
% nH = 1;
% K_H = 1.5*1e-3; %1.5 microM
V_NaCa = (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out/Ca_in))*1e3; %mV
I_NaCa_gall_susa = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa);
% gall_susa_fit = lsqcurvefit(I_NaCa_gall_susa,-0.07,Vm_data,I_data);
plot(Vm,I_NaCa_gall_susa(gall_susa_fit,Vm),'--','LineWidth',4);

%testing flanagan et al 2018 model
I_NaCa_flanagan = @(x,Vm) x*((Na_in.^3)./(Na_out.^3).*(exp(0.5*Vm*1e-3*F/(R*T))) ...
   - (Ca_in)./(Ca_out).*(exp(0.5*Vm*1e-3*F/(R*T)))) ;
% flanagan_fit = lsqcurvefit(I_NaCa_flanagan,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_flanagan(flanagan_fit,Vm),':','LineWidth',4);

%testing my model
nH2 = 1.6;
kH2 = 20;
I_NaCa_my1 = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH) ...
    .*(Na_in.^nH2)./(Na_in.^nH2 + kH2^nH2).*(Vm-V_NaCa);
% my_fit = lsqcurvefit(I_NaCa_my,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_my1(Cainhill_and_Nainhill,Vm),':','LineWidth',4);

I_NaCa_linear = @(x,Vm) x.*(Vm-V_NaCa);
% my_fit2 = lsqcurvefit(I_NaCa_my2,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_linear(linear_fit,Vm),':','LineWidth',4);


lgd = legend('Matsuoka & Hilgemann 1992 Fig9A data','Egger & Niggli 1999 model fit',...
    'Kimura et al. 1987 model fit', 'Gall & Susa 1999 model fit',...
    'Flanagan et al. 2018 model fit', 'Ca_i Hill + Na_i Hill', 'Linear');
lgd.Location = 'northwest';

line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-50,50],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

title('[Na^+]_i = 50 mM, [Ca^{2+}]_i = 1 \muM');
xlim([-150,100]); ylim([-50,100])


%data from Fig9A - 20 mM Nain
[V_data, I_data] = csvimport('data/matsuoka_hilgemann_1992_fig9a_20.csv',...
                    'columns', [1, 2] ,'noHeader', true);              
Na_in = 20;%mM          
figure(3); clf;
plot(Vm_data,I_data,'o','LineWidth',4,'MarkerSize',10); hold on
xlabel('V_m (mV)'); ylabel('current, pA'); set(gca,'FontSize',20);

%testing egger & niggli 1999 model
%x = [0.2; 0.5;]; %A l^4 /mol^4, partition coefficient (symmetric reversal potential)
I_NaCa_egger = @(x,Vm) x*((Na_in.*1e-3).^3.*(Ca_out.*1e-3).*(exp(0.5*Vm*1e-3*F/(R*T))) ...
   - ((Na_out.*1e-3).^3.*(Ca_in.*1e-3).*(exp(-0.5*Vm*1e-3*F/(R*T)))));
% egger_fit = lsqcurvefit(I_NaCa_egger,2e3,Vm_data,I_data);
plot(Vm,I_NaCa_egger(egger_fit,Vm),'-','LineWidth',4);


%testing kimura et al 1987 model
% x = [2.07e-5; 0.38;]; %muA muF /mM^4, partition coefficient (symmetric reversal potential)
I_NaCa_kimura = @(x,Vm) x*((Na_in.^3).*(Ca_out).*(exp(0.38*Vm*1e-3*F/(R*T))) ...
   - (Na_out.^3).*(Ca_in).*(exp(-(1-0.38)*Vm*1e-3*F/(R*T)))) ;
% kimura_fit = lsqcurvefit(I_NaCa_kimura,2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_kimura(kimura_fit,Vm),'-.','LineWidth',4);

%testing gall & susa 1999 model
% nH = 1;
% K_H = 1.5*1e-3; %1.5 microM
V_NaCa = (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out/Ca_in))*1e3; %mV
I_NaCa_gall_susa = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa);
% gall_susa_fit = lsqcurvefit(I_NaCa_gall_susa,-0.07,Vm_data,I_data);
plot(Vm,I_NaCa_gall_susa(gall_susa_fit,Vm),'--','LineWidth',4);

%testing flanagan et al 2018 model
I_NaCa_flanagan = @(x,Vm) x*((Na_in.^3)./(Na_out.^3).*(exp(0.5*Vm*1e-3*F/(R*T))) ...
   - (Ca_in)./(Ca_out).*(exp(0.5*Vm*1e-3*F/(R*T)))) ;
% flanagan_fit = lsqcurvefit(I_NaCa_flanagan,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_flanagan(flanagan_fit,Vm),':','LineWidth',4);


%testing my model
nH2 = 1.6;
kH2 = 20;
I_NaCa_my1 = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH) ...
    .*(Na_in.^nH2)./(Na_in.^nH2 + kH2^nH2).*(Vm-V_NaCa);
% my_fit = lsqcurvefit(I_NaCa_my,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_my1(Cainhill_and_Nainhill,Vm),':','LineWidth',4);

I_NaCa_linear = @(x,Vm) x.*(Vm-V_NaCa);
% my_fit2 = lsqcurvefit(I_NaCa_my2,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_linear(linear_fit,Vm),':','LineWidth',4);


lgd = legend('Matsuoka & Hilgemann 1992 Fig9A data','Egger & Niggli 1999 model fit',...
    'Kimura et al. 1987 model fit', 'Gall & Susa 1999 model fit',...
    'Flanagan et al. 2018 model fit', 'Ca_i Hill + Na_i Hill', 'Linear');
lgd.Location = 'northwest';

line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-50,50],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

title('[Na^+]_i = 20 mM, [Ca^{2+}]_i = 1 \muM');
xlim([-150,100]); ylim([-50,100])

%data from Fig9A - 5 mM Nain
[V_data, I_data] = csvimport('data/matsuoka_hilgemann_1992_fig9a_5.csv',...
                    'columns', [1, 2] ,'noHeader', true);              
Na_in = 5;%mM          
figure(4); clf;
plot(Vm_data,I_data,'o','LineWidth',4,'MarkerSize',10); hold on
xlabel('V_m (mV)'); ylabel('current, pA'); set(gca,'FontSize',20);

%testing egger & niggli 1999 model
%x = [0.2; 0.5;]; %A l^4 /mol^4, partition coefficient (symmetric reversal potential)
I_NaCa_egger = @(x,Vm) x*((Na_in.*1e-3).^3.*(Ca_out.*1e-3).*(exp(0.5*Vm*1e-3*F/(R*T))) ...
   - ((Na_out.*1e-3).^3.*(Ca_in.*1e-3).*(exp(-0.5*Vm*1e-3*F/(R*T)))));
% egger_fit = lsqcurvefit(I_NaCa_egger,2e3,Vm_data,I_data);
plot(Vm,I_NaCa_egger(egger_fit,Vm),'-','LineWidth',4);


%testing kimura et al 1987 model
% x = [2.07e-5; 0.38;]; %muA muF /mM^4, partition coefficient (symmetric reversal potential)
I_NaCa_kimura = @(x,Vm) x*((Na_in.^3).*(Ca_out).*(exp(0.38*Vm*1e-3*F/(R*T))) ...
   - (Na_out.^3).*(Ca_in).*(exp(-(1-0.38)*Vm*1e-3*F/(R*T)))) ;
% kimura_fit = lsqcurvefit(I_NaCa_kimura,2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_kimura(kimura_fit,Vm),'-.','LineWidth',4);

%testing gall & susa 1999 model
% nH = 1;
% K_H = 1.5*1e-3; %1.5 microM
V_NaCa = (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out/Ca_in))*1e3; %mV
I_NaCa_gall_susa = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa);
% gall_susa_fit = lsqcurvefit(I_NaCa_gall_susa,-0.07,Vm_data,I_data);
plot(Vm,I_NaCa_gall_susa(gall_susa_fit,Vm),'--','LineWidth',4);

%testing flanagan et al 2018 model
I_NaCa_flanagan = @(x,Vm) x*((Na_in.^3)./(Na_out.^3).*(exp(0.5*Vm*1e-3*F/(R*T))) ...
   - (Ca_in)./(Ca_out).*(exp(0.5*Vm*1e-3*F/(R*T)))) ;
% flanagan_fit = lsqcurvefit(I_NaCa_flanagan,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_flanagan(flanagan_fit,Vm),':','LineWidth',4);

%testing my model
nH2 = 1.6;
kH2 = 20;
I_NaCa_my1 = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH) ...
    .*(Na_in.^nH2)./(Na_in.^nH2 + kH2^nH2).*(Vm-V_NaCa);
% my_fit = lsqcurvefit(I_NaCa_my,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_my1(Cainhill_and_Nainhill,Vm),':','LineWidth',4);

I_NaCa_linear = @(x,Vm) x.*(Vm-V_NaCa);
% my_fit2 = lsqcurvefit(I_NaCa_my2,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_linear(linear_fit,Vm),':','LineWidth',4);

lgd = legend('Matsuoka & Hilgemann 1992 Fig9A data','Egger & Niggli 1999 model fit',...
    'Kimura et al. 1987 model fit', 'Gall & Susa 1999 model fit',...
    'Flanagan et al. 2018 model fit', 'Ca_i Hill + Na_i Hill', 'Linear');
lgd.Location = 'northwest';

line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-50,50],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

title('[Na^+]_i = 5 mM, [Ca^{2+}]_i = 1 \muM');
xlim([-150,100]); ylim([-50,100])

%Fig 9B - I-V curves
%recording at 34-36 deg C (described as exactly the same method as
%Hilgemann et al. 1992)
%recording of an outward exchange current at various membrane potentials
%and various internal sodium concentrations
%Na_in = 100; %mM varied here
Ca_in = 0.1e-3; %mM (1 microM)
Na_out = 150; %mM
Ca_out = 2; %mM
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 308.15; %K, 35 deg C

%Vm for I-V curve
Vm = (-150:5:100); %mV

%data from Fig9b - first points
[Vm_data, I_data] = csvimport('data/matsuoka_hilgemann_1992_fig9b_100.csv',...
                    'columns', [1, 2] ,'noHeader', true);
Na_in = 100;%mM                
figure(5); clf;
plot(Vm_data,I_data,'o','LineWidth',4,'MarkerSize',10); hold on
xlabel('V_m (mV)'); ylabel('current, pA'); set(gca,'FontSize',20);

%testing egger & niggli 1999 model
%x = [0.2; 0.5;]; %A l^4 /mol^4, partition coefficient (symmetric reversal potential)
I_NaCa_egger = @(x,Vm) x*((Na_in.*1e-3).^3.*(Ca_out.*1e-3).*(exp(0.5*Vm*1e-3*F/(R*T))) ...
   - ((Na_out.*1e-3).^3.*(Ca_in.*1e-3).*(exp(-0.5*Vm*1e-3*F/(R*T)))));
% egger_fit = lsqcurvefit(I_NaCa_egger,2e4,Vm_data,I_data);
plot(Vm,I_NaCa_egger(egger_fit,Vm),'-','LineWidth',4);


%testing kimura et al 1987 model
% x = [2.07e-5; 0.38;]; %muA muF /mM^4, partition coefficient (symmetric reversal potential)
I_NaCa_kimura = @(x,Vm) x*((Na_in.^3).*(Ca_out).*(exp(0.38*Vm*1e-3*F/(R*T))) ...
   - (Na_out.^3).*(Ca_in).*(exp(-(1-0.38)*Vm*1e-3*F/(R*T)))) ;
% kimura_fit = lsqcurvefit(I_NaCa_kimura,2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_kimura(kimura_fit,Vm),'-.','LineWidth',4);

%testing gall & susa 1999 model
%hill params for gall & susa
% nH = 1;
% K_H = 1.5*1e-3; %1.5 microM
V_NaCa = (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out/Ca_in))*1e3; %mV
I_NaCa_gall_susa = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa);
% gall_susa_fit = lsqcurvefit(I_NaCa_gall_susa,-0.07,Vm_data,I_data);
plot(Vm,I_NaCa_gall_susa(gall_susa_fit,Vm),'--','LineWidth',4);

%testing flanagan et al 2018 model
I_NaCa_flanagan = @(x,Vm) x*((Na_in.^3)./(Na_out.^3).*(exp(0.5*Vm*1e-3*F/(R*T))) ...
   - (Ca_in)./(Ca_out).*(exp(0.5*Vm*1e-3*F/(R*T)))) ;
% flanagan_fit = lsqcurvefit(I_NaCa_flanagan,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_flanagan(flanagan_fit,Vm),':','LineWidth',4);

%testing my model
nH2 = 1.6;
kH2 = 20;
I_NaCa_my1 = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH) ...
    .*(Na_in.^nH2)./(Na_in.^nH2 + kH2^nH2).*(Vm-V_NaCa);
% my_fit = lsqcurvefit(I_NaCa_my,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_my1(Cainhill_and_Nainhill,Vm),':','LineWidth',4);

I_NaCa_linear = @(x,Vm) x.*(Vm-V_NaCa);
% my_fit2 = lsqcurvefit(I_NaCa_my2,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_linear(linear_fit,Vm),':','LineWidth',4);

lgd = legend('Matsuoka & Hilgemann 1992 Fig9B data','Egger & Niggli 1999 model fit',...
    'Kimura et al. 1987 model fit', 'Gall & Susa 1999 model fit',...
    'Flanagan et al. 2018 model fit', 'Ca_i Hill + Na_i Hill', 'Linear');
lgd.Location = 'northwest';


line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-50,50],'LineStyle','--', 'Color', 'k','HandleVisibility','off');


title('[Na^+]_i = 100 mM, [Ca^{2+}]_i = 0.1 \muM');
xlim([-150,100]); ylim([-50,100])


%data from Fig9b -more points
[Vm_data, I_data] = csvimport('data/matsuoka_hilgemann_1992_fig9b_10.csv',...
                    'columns', [1, 2] ,'noHeader', true);
Na_in = 10;%mM                
figure(6); clf;
plot(Vm_data,I_data,'o','LineWidth',4,'MarkerSize',10); hold on
xlabel('V_m (mV)'); ylabel('current, pA'); set(gca,'FontSize',20);

%testing egger & niggli 1999 model
%x = [0.2; 0.5;]; %A l^4 /mol^4, partition coefficient (symmetric reversal potential)
I_NaCa_egger = @(x,Vm) x*((Na_in.*1e-3).^3.*(Ca_out.*1e-3).*(exp(0.5*Vm*1e-3*F/(R*T))) ...
   - ((Na_out.*1e-3).^3.*(Ca_in.*1e-3).*(exp(-0.5*Vm*1e-3*F/(R*T)))));
% egger_fit = lsqcurvefit(I_NaCa_egger,2e4,Vm_data,I_data);
plot(Vm,I_NaCa_egger(egger_fit,Vm),'-','LineWidth',4);


%testing kimura et al 1987 model
% x = [2.07e-5; 0.38;]; %muA muF /mM^4, partition coefficient (symmetric reversal potential)
I_NaCa_kimura = @(x,Vm) x*((Na_in.^3).*(Ca_out).*(exp(0.38*Vm*1e-3*F/(R*T))) ...
   - (Na_out.^3).*(Ca_in).*(exp(-(1-0.38)*Vm*1e-3*F/(R*T)))) ;
% kimura_fit = lsqcurvefit(I_NaCa_kimura,2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_kimura(kimura_fit,Vm),'-.','LineWidth',4);

%testing gall & susa 1999 model
%hill params for gall & susa
% nH = 1;
% K_H = 1.5*1e-3; %1.5 microM
V_NaCa = (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out/Ca_in))*1e3; %mV
I_NaCa_gall_susa = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa);
% gall_susa_fit = lsqcurvefit(I_NaCa_gall_susa,-0.07,Vm_data,I_data);
plot(Vm,I_NaCa_gall_susa(gall_susa_fit,Vm),'--','LineWidth',4);

%testing flanagan et al 2018 model
I_NaCa_flanagan = @(x,Vm) x*((Na_in.^3)./(Na_out.^3).*(exp(0.5*Vm*1e-3*F/(R*T))) ...
   - (Ca_in)./(Ca_out).*(exp(0.5*Vm*1e-3*F/(R*T)))) ;
% flanagan_fit = lsqcurvefit(I_NaCa_flanagan,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_flanagan(flanagan_fit,Vm),':','LineWidth',4);

%testing my model
nH2 = 1.6;
kH2 = 20;
I_NaCa_my1 = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH) ...
    .*(Na_in.^nH2)./(Na_in.^nH2 + kH2^nH2).*(Vm-V_NaCa);
% my_fit = lsqcurvefit(I_NaCa_my,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_my1(Cainhill_and_Nainhill,Vm),':','LineWidth',4);

I_NaCa_linear = @(x,Vm) x.*(Vm-V_NaCa);
% my_fit2 = lsqcurvefit(I_NaCa_my2,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_linear(linear_fit,Vm),':','LineWidth',4);

lgd = legend('Matsuoka & Hilgemann 1992 Fig9B data','Egger & Niggli 1999 model fit',...
    'Kimura et al. 1987 model fit', 'Gall & Susa 1999 model fit',...
    'Flanagan et al. 2018 model fit', 'Ca_i Hill + Na_i Hill', 'Linear');
lgd.Location = 'northwest';


line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-50,50],'LineStyle','--', 'Color', 'k','HandleVisibility','off');


title('[Na^+]_i = 10 mM, [Ca^{2+}]_i = 0.1 \muM');
xlim([-150,100]); ylim([-50,100])

