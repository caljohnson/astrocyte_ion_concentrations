%Hilgemann et al. 1992 data test
%testing each model against the data and I-V curves

%Fig 2C - I-V curves
%recording at 34-36 deg C
%recording of an outward exchange current at various membrane potentials
%pipette solution (external): 0mM Na, 8mM Ca
%cytoplasmic solution (internal): 100mM Na, 2 microM Ca
Na_in = 100; %mM
Ca_in = 2e-3; %mM
Na_out = 0; %mM
Ca_out = 8; %mM
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 308.15; %K, 35 deg C

%Vm for I-V curve
Vm = (-120:5:60); %mV

%data from Fig2C
[Vm_data, I_data] = csvimport('data/hilgemann_1992_fig2c.csv',...
                    'columns', [1, 2] ,'noHeader', true);
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
plot(Vm,I_NaCa_kimura(kimura_fit,Vm),'-','LineWidth',4);


%testing flanagan et al 2018 model
I_NaCa_flanagan = @(x,Vm) x*((Na_in.^3)./(Na_out.^3).*(exp(0.5*Vm*1e-3*F/(R*T))) ...
   - (Ca_in)./(Ca_out).*(exp(0.5*Vm*1e-3*F/(R*T)))) ;
%Flanagan model is undefined/infinite for Na_out = 0 mM
% flanagan_fit = lsqcurvefit(I_NaCa_flanagan,-2e-5,Vm_data,I_data);
% plot(Vm,I_NaCa_flanagan(flanagan_fit,Vm),'-','LineWidth',4);

%testing gall & susa 1999 model
nH = 1;
K_H = 1.5*1e-3; %1.5 microM
V_NaCa = (R*T/F)*(3*log(Na_out/Na_in)-log(Ca_out/Ca_in))*1e3; %mV
V_NaCa = -150; %since Na_out = 0, V_NaCa = -inf
I_NaCa_gall_susa = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa);
gall_susa_fit = lsqcurvefit(I_NaCa_gall_susa,-0.07,Vm_data,I_data);
plot(Vm,I_NaCa_gall_susa(gall_susa_fit,Vm),'-','LineWidth',4);

lgd = legend('Hilgemann et al. 1992 Fig2C data','Egger & Niggli 1999 model fit',...
    'Kimura et al. 1987 model fit', 'Gall & Susa 1999 model fit');
lgd.Location = 'northwest';

line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[0,25],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

%data from Hilgemann et al. 1992, Fig7B
[V_data2, I_data2] = csvimport('data/hilgemann_1992_fig7b.csv',...
                    'columns', [1, 2] ,'noHeader', true);              
          
figure(2); clf;
plot(V_data2,I_data2,'o','LineWidth',4,'MarkerSize',10); hold on
xlabel('V_m (mV)'); ylabel('current, pA'); set(gca,'FontSize',20);

%testing egger & niggli 1999 model
%x = [0.2; 0.5;]; %A l^4 /mol^4, partition coefficient (symmetric reversal potential)
I_NaCa_egger2 = @(x,Vm) x*((Na_in.*1e-3).^3.*(Ca_out.*1e-3).*(exp(0.5*Vm*1e-3*F./(R*T))) ...
   - ((Na_out.*1e-3).^3.*(Ca_in.*1e-3).*(exp(-0.5*Vm*1e-3*F./(R*T)))));
egger_fit2 = lsqcurvefit(I_NaCa_egger2,2e4,V_data2,I_data2);
plot(Vm,I_NaCa_egger2(egger_fit2,Vm),'-','LineWidth',4);


%testing kimura et al 1987 model
% x = [2.07e-5; 0.38;]; %muA muF /mM^4, partition coefficient (symmetric reversal potential)
I_NaCa_kimura2 = @(x,Vm) x*((Na_in.^3).*(Ca_out).*(exp(0.38*Vm*1e-3*F./(R*T))) ...
   - (Na_out.^3).*(Ca_in).*(exp(-(1-0.38)*Vm*1e-3*F./(R*T)))) ;
kimura_fit2 = lsqcurvefit(I_NaCa_kimura2,2e-5,V_data2,I_data2);
plot(Vm,I_NaCa_kimura2(kimura_fit2,Vm),'-','LineWidth',4);


%testing flanagan et al 2018 model
I_NaCa_flanagan2 = @(x,Vm) x*((Na_in.^3)./(Na_out.^3).*(exp(0.5*Vm*1e-3*F./(R*T))) ...
   - (Ca_in)./(Ca_out).*(exp(0.5*Vm*1e-3*F./(R*T)))) ;
%Flanagan model is undefined/infinite for Na_out = 0 mM
% flanagan_fit = lsqcurvefit(I_NaCa_flanagan,-2e-5,Vm_data,I_data);
% plot(Vm,I_NaCa_flanagan(flanagan_fit,Vm),'-','LineWidth',4);

% testing gall & susa 1999 model
nH = 1;
K_H = 1.5*1e-3; %1.5 microM
V_NaCa2 = (R*T/F).*(3*log(Na_out/Na_in)-log(Ca_out/Ca_in))*1e3; %mV
V_NaCa = -150; %since Na_out = 0, V_NaCa = -inf
I_NaCa_gall_susa2 = @(x,Vm) x*(Ca_in.^nH)./(Ca_in.^nH + K_H^nH).*(Vm-V_NaCa);
gall_susa_fit2 = lsqcurvefit(I_NaCa_gall_susa,-0.07,V_data2,I_data2);
plot(Vm,I_NaCa_gall_susa2(gall_susa_fit2,Vm),'-','LineWidth',4);

lgd = legend('Hilgemann et al. 1992 Fig7B data','Egger & Niggli 1999 model fit',...
    'Kimura et al. 1987 model fit', 'Gall & Susa 1999 model fit');
lgd.Location = 'northwest';

line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
line([0,0],[0,60],'LineStyle','--', 'Color', 'k','HandleVisibility','off');

