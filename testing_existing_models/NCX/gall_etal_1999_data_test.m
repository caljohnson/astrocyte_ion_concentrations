%Gall et al. 1999 data test
%testing each model against the data and I-V curves

%Fig 4A - I-V curves
%recording at ? deg C 
%recording of exchange current at various membrane potentials
%and various external calcium concentrations
%Ca_out = ; %mM varied here
Ca_in = 14e-3; %mM (14 microM)
Na_out = 138; %mM
Na_in = 30; %mM
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 308.15; %K, 35 deg C

%Vm for I-V curve
Vm = (-80:1:40); %mV

%data from Fig4A - first points
[Vm_data, I_data] = csvimport('data/gall_etal_fig4A_Cae_5mM.csv',...
                    'columns', [1, 2] ,'noHeader', true);
Ca_out = 5;%mM                
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
% nH = 5;
% K_H = 1.5*1e-3; %1.5 microM
nH = 1;
% K_H = 0.05*1e-3;%mM
V_NaCa = (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out/Ca_in))*1e3; %mV
I_NaCa_gall_susa = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa);
gall_susa_fit = lsqcurvefit(I_NaCa_gall_susa,-0.07,Vm_data,I_data);
plot(Vm,I_NaCa_gall_susa(gall_susa_fit,Vm),'--','LineWidth',4);

%testing flanagan et al 2018 model
I_NaCa_flanagan = @(x,Vm) x*((Na_in.^3)./(Na_out.^3).*(exp(0.5*Vm*1e-3*F/(R*T))) ...
   - (Ca_in)./(Ca_out).*(exp(0.5*Vm*1e-3*F/(R*T)))) ;
flanagan_fit = lsqcurvefit(I_NaCa_flanagan,-2e-5,Vm_data,I_data);
plot(Vm,I_NaCa_flanagan(flanagan_fit,Vm),':','LineWidth',4);

% %testing my model
% nH2 = 1.6;
% kH2 = 20;
% I_NaCa_my1 = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH) ...
%     .*(Na_in.^nH2)./(Na_in.^nH2 + kH2^nH2).*(Vm-V_NaCa);
% Cainhill_and_Nainhill = lsqcurvefit(I_NaCa_my1,-2e-5,Vm_data,I_data);
% plot(Vm,I_NaCa_my1(Cainhill_and_Nainhill,Vm),':','LineWidth',4);
% 
% I_NaCa_linear = @(x,Vm) x.*(Vm-V_NaCa);
% linear_fit = lsqcurvefit(I_NaCa_linear,2,Vm_data,I_data);
% plot(Vm,I_NaCa_linear(linear_fit,Vm),':','LineWidth',4);


lgd = legend('Gall et al. 1999 Fig4A data','Egger & Niggli 1999 model fit',...
    'Kimura et al. 1987 model fit', 'Gall & Susa 1999 model fit',...
    'Flanagan et al. 2018 model fit');%, 'Ca_i Hill + Na_i Hill', 'Linear');
lgd.Location = 'southeast';

line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-50,50],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

title('[Ca^{2+}]_e = 5 mM');
xlim([-80,40]); ylim([-5,5])


%data from Fig4A - 2.6 mM Caout
[Vm_data, I_data] = csvimport('data/gall_etal_fig4A_Cae_2p6mM.csv',...
                    'columns', [1, 2] ,'noHeader', true);              
Ca_out = 2.6;%mM          
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

% %testing my model
% nH2 = 1.6;
% kH2 = 20;
% I_NaCa_my1 = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH) ...
%     .*(Na_in.^nH2)./(Na_in.^nH2 + kH2^nH2).*(Vm-V_NaCa);
% % my_fit = lsqcurvefit(I_NaCa_my,-2e-5,Vm_data,I_data);
% plot(Vm,I_NaCa_my1(Cainhill_and_Nainhill,Vm),':','LineWidth',4);
% 
% I_NaCa_linear = @(x,Vm) x.*(Vm-V_NaCa);
% % my_fit2 = lsqcurvefit(I_NaCa_my2,-2e-5,Vm_data,I_data);
% plot(Vm,I_NaCa_linear(linear_fit,Vm),':','LineWidth',4);


lgd = legend('Gall et al. 1999 Fig4A data','Egger & Niggli 1999 model fit',...
    'Kimura et al. 1987 model fit', 'Gall & Susa 1999 model fit',...
    'Flanagan et al. 2018 model fit');%, 'Ca_i Hill + Na_i Hill', 'Linear');
lgd.Location = 'southeast';

line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-50,50],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

title('[Ca^{2+}]_e = 2.6 mM');
xlim([-80,40]); ylim([-5,5])


%data from Fig9A - 1.6mM  Caout
[Vm_data, I_data] = csvimport('data/gall_etal_fig4A_Cae_1p6mM.csv',...
                    'columns', [1, 2] ,'noHeader', true);              
Ca_out = 1.6;%mM          
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


% %testing my model
% nH2 = 1.6;
% kH2 = 20;
% I_NaCa_my1 = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH) ...
%     .*(Na_in.^nH2)./(Na_in.^nH2 + kH2^nH2).*(Vm-V_NaCa);
% % my_fit = lsqcurvefit(I_NaCa_my,-2e-5,Vm_data,I_data);
% plot(Vm,I_NaCa_my1(Cainhill_and_Nainhill,Vm),':','LineWidth',4);
% 
% I_NaCa_linear = @(x,Vm) x.*(Vm-V_NaCa);
% % my_fit2 = lsqcurvefit(I_NaCa_my2,-2e-5,Vm_data,I_data);
% plot(Vm,I_NaCa_linear(linear_fit,Vm),':','LineWidth',4);
% 

lgd = legend('Gall et al. 1999 Fig4A data','Egger & Niggli 1999 model fit',...
    'Kimura et al. 1987 model fit', 'Gall & Susa 1999 model fit',...
    'Flanagan et al. 2018 model fit');%, 'Ca_i Hill + Na_i Hill', 'Linear');
lgd.Location = 'southeast';

line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-50,50],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

title('[Ca^{2+}]_e = 1.6 mM');
xlim([-80,40]); ylim([-5,4])

%data from Fig4A - 1 mM Caout
[Vm_data, I_data] = csvimport('data/gall_etal_fig4A_Cae_1mM.csv',...
                    'columns', [1, 2] ,'noHeader', true);              
Ca_out = 1;%mM          
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

% %testing my model
% nH2 = 1.6;
% kH2 = 20;
% I_NaCa_my1 = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH) ...
%     .*(Na_in.^nH2)./(Na_in.^nH2 + kH2^nH2).*(Vm-V_NaCa);
% % my_fit = lsqcurvefit(I_NaCa_my,-2e-5,Vm_data,I_data);
% plot(Vm,I_NaCa_my1(Cainhill_and_Nainhill,Vm),':','LineWidth',4);
% 
% I_NaCa_linear = @(x,Vm) x.*(Vm-V_NaCa);
% % my_fit2 = lsqcurvefit(I_NaCa_my2,-2e-5,Vm_data,I_data);
% plot(Vm,I_NaCa_linear(linear_fit,Vm),':','LineWidth',4);

lgd = legend('Gall et al. 1999 Fig4A data','Egger & Niggli 1999 model fit',...
    'Kimura et al. 1987 model fit', 'Gall & Susa 1999 model fit',...
    'Flanagan et al. 2018 model fit');%, 'Ca_i Hill + Na_i Hill', 'Linear');
lgd.Location = 'northwest';

line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-50,50],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

title('[Ca^{2+}]_e = 1 mM');
xlim([-80,40]); ylim([-5,5])

%Fig 4B - I-V curves
%recording at ? deg C (described as exactly the same method as
%Hilgemann et al. 1992)

%data from Fig4b - first points
[Vm_data, I_data] = csvimport('data/gall_etal_fig4B_Cai_14muM.csv',...
                    'columns', [1, 2] ,'noHeader', true);
Ca_out = 2.6; %MM
Ca_in = 14e-3;%mM (14 microM)                
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

% %testing my model
% nH2 = 1.6;
% kH2 = 20;
% I_NaCa_my1 = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH) ...
%     .*(Na_in.^nH2)./(Na_in.^nH2 + kH2^nH2).*(Vm-V_NaCa);
% % my_fit = lsqcurvefit(I_NaCa_my,-2e-5,Vm_data,I_data);
% plot(Vm,I_NaCa_my1(Cainhill_and_Nainhill,Vm),':','LineWidth',4);
% 
% I_NaCa_linear = @(x,Vm) x.*(Vm-V_NaCa);
% % my_fit2 = lsqcurvefit(I_NaCa_my2,-2e-5,Vm_data,I_data);
% plot(Vm,I_NaCa_linear(linear_fit,Vm),':','LineWidth',4);

lgd = legend('Gall et al. 1999 Fig4B data','Egger & Niggli 1999 model fit',...
    'Kimura et al. 1987 model fit', 'Gall & Susa 1999 model fit',...
    'Flanagan et al. 2018 model fit');%, 'Ca_i Hill + Na_i Hill', 'Linear');
lgd.Location = 'southeast';


line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-50,50],'LineStyle','--', 'Color', 'k','HandleVisibility','off');


title('[Ca^{2+}]_i = 14 \muM');
xlim([-80,40]); ylim([-5,5])


%data from Fig4b -more points
[Vm_data, I_data] = csvimport('data/gall_etal_fig4B_Cai_3p5muM.csv',...
                    'columns', [1, 2] ,'noHeader', true);
Ca_in = 3.5e-3;%mM (3.5 microM)                 
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

% %testing my model
% nH2 = 1.6;
% kH2 = 20;
% I_NaCa_my1 = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH) ...
%     .*(Na_in.^nH2)./(Na_in.^nH2 + kH2^nH2).*(Vm-V_NaCa);
% % my_fit = lsqcurvefit(I_NaCa_my,-2e-5,Vm_data,I_data);
% plot(Vm,I_NaCa_my1(Cainhill_and_Nainhill,Vm),':','LineWidth',4);
% 
% I_NaCa_linear = @(x,Vm) x.*(Vm-V_NaCa);
% % my_fit2 = lsqcurvefit(I_NaCa_my2,-2e-5,Vm_data,I_data);
% plot(Vm,I_NaCa_linear(linear_fit,Vm),':','LineWidth',4);

lgd = legend('Gall et al. 1999 Fig4B data','Egger & Niggli 1999 model fit',...
    'Kimura et al. 1987 model fit', 'Gall & Susa 1999 model fit',...
    'Flanagan et al. 2018 model fit');%, 'Ca_i Hill + Na_i Hill', 'Linear');
lgd.Location = 'southeast';


line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[-50,50],'LineStyle','--', 'Color', 'k','HandleVisibility','off');


title('[Ca^{2+}]_i = 3.5 \muM');
xlim([-80,40]); ylim([-5,5])

