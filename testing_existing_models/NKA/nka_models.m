%NKA models

%Huguet et al. 2016
p=1;
kK = 2; %mM, half activation [K+]e concentration for the astrocytic NKA
kNa = 7.7; %mM, half activation [Na+]i concentration for the astrocytic NKA
I_NKA_huguet = @(K_out, Na_in) p.*(K_out./(kK + K_out)).^2.*(Na_in./(kNa + Na_in)).^3;

%Kager et al. 2000
kK = 3.5; %mM, half activation [K+]e concentration for the astrocytic NKA
kNa = 10; %mM, half activation [Na+]i concentration for the astrocytic NKA
I_NKA_kager = @(K_out, Na_in) (K_out./(kK + K_out)).^2.*(Na_in./(kNa + Na_in)).^3;


%Cressman et al. 2009
kNa = 25.0;
kK  = 5.5;
I_NKA_cressman = @(K_out, Na_in) p.*(1+exp((kNa-Na_in)/3)).^(-1)...
                        .*(1+exp(kK-K_out)).^(-1);

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
plot(K_outs,I_NKA_huguet(K_outs,Na_in),'b-','LineWidth',4); hold on;
plot(K_outs,I_NKA_kager(K_outs,Na_in),'r--','LineWidth',4);
plot(K_outs,I_NKA_cressman(K_outs,Na_in),'m:','LineWidth',4); 
% line([-150,100],[0,0],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
% line([0,0],[-2,2],'LineStyle','--', 'Color', 'k','HandleVisibility','off');
xlabel('[K^+]_e (mM)'); ylabel('current');set(gca,'FontSize',20);
legend('Huguet et al. 2016 model', 'Kager et al. 2000 model', 'Cressman et al. 2009 model');

Na_ins = linspace(0,100);
figure(2);clf
plot(Na_ins,I_NKA_huguet(K_out,Na_ins),'b-','LineWidth',4); hold on
plot(Na_ins,I_NKA_kager(K_out,Na_ins),'r--','LineWidth',4);
plot(Na_ins,I_NKA_cressman(K_out,Na_ins),'m:','LineWidth',4);
xlabel('[Na^+]_i (mM)'); ylabel('current');set(gca,'FontSize',20);
legend('Huguet et al. 2016 model', 'Kager et al. 2000 model', 'Cressman et al. 2009 model');


%vary kK
K_outs = linspace(0,100);
figure(3);clf
kKs = linspace(1,21,10);
lbls = {};
for jj=1:10
    kK = kKs(jj);
%     I_NKA = @(K_out, Na_in) p.*(K_out./(kK + K_out)).^2.*(Na_in./(kNa + Na_in)).^3;
    I_NKA = @(K_out, Na_in) p.*(1+exp((kNa-Na_in)/3)).^(-1)...
                        .*(1+exp(kK-K_out)).^(-1);
    plot(K_outs, I_NKA(K_outs,Na_in),'LineWidth',4); hold on;
    lbls{jj} = ['k_K = ' num2str(kK) ' mM'];
end
xlabel('[K^+]_e (mM)'); ylabel('current');set(gca,'FontSize',20);
legend(lbls);


%vary kNa
Na_ins = linspace(0,100);
figure(4);clf
kNas = linspace(1,100,10);
for jj=1:10
    kNa = kNas(jj);
    I_NKA = @(K_out, Na_in) p.*(1+exp((kNa-Na_in)/3)).^(-1)...
                        .*(1+exp(kK-K_out)).^(-1);
    plot(Na_ins, I_NKA(K_out,Na_ins),'LineWidth',4); hold on;
    lbls2{jj} = ['k_{Na} = ' num2str(kNa) ' mM'];
end
xlabel('[Na^+]_i (mM)'); ylabel('current');set(gca,'FontSize',20);
legend(lbls2);