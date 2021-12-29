%Nakao & Gadsby 1989 data test
%of NKA models

%constants
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 298.15; %K, absolute temperature (25 C in Kirischuk et al. 2012)

%data from Figure 8
[Na_in_data, I_data] = csvimport('nakao_gadsby_1989_fig8.csv',...
                    'columns', [1, 2] ,'noHeader', true);
%renormalize Idata since its weird
I_data = I_data./max(abs(I_data));
% Na_out = 150; %mM
K_out = 5.4; %mM
% V = 0; %mV

figure(1); clf;
plot(Na_in_data, I_data, 'o', 'MarkerSize',10, 'LineWidth',4); hold on

%Huguet et al. 2016
p0=300;
kK = 2; %mM, half activation [K+]e concentration for the astrocytic NKA
kNa = 7.7; %mM, half activation [Na+]i concentration for the astrocytic NKA
% I_NKA_huguet = @(K_out, Na_in) p0.*(K_out./(kK + K_out)).^2.*(Na_in./(kNa + Na_in)).^3;
I_NKA_huguet = @(K_out, Na_in) (K_out./(kK + K_out)).^2.*(Na_in./(kNa + Na_in)).^3;

%Kager et al. 2000
p1 = 500;
kK = 3.5; %mM, half activation [K+]e concentration for the astrocytic NKA
kNa = 10; %mM, half activation [Na+]i concentration for the astrocytic NKA
% I_NKA_kager = @(K_out, Na_in) p1.*(K_out./(kK + K_out)).^2.*(Na_in./(kNa + Na_in)).^3;
I_NKA_kager = @(K_out, Na_in) (K_out./(kK + K_out)).^2.*(Na_in./(kNa + Na_in)).^3;


%Cressman et al. 2009
p2 = 240;
kNa = 25.0;
kK  = 5.5;
% I_NKA_cressman = @(K_out, Na_in) p2.*(1+exp((kNa-Na_in)/3)).^(-1)...
%                         .*(1+exp(kK-K_out)).^(-1);
I_NKA_cressman = @(K_out, Na_in) (1+exp((kNa-Na_in)/3)).^(-1)...
                        .*(1+exp(kK-K_out)).^(-1);
                    
%Oschmann et al. 2017
kNa = 10; %mM
kK  = 1.5; %mM
I_NKA_oschmann = @(K_out, Na_in) ((Na_in.^(1.5))./(Na_in.^(1.5) + kNa^(1.5)))...
    .*(K_out./(K_out + kK));

%my model
I_NKA_Nahill = @(x, Na_in)  (Na_in.^x(2)./(x(1).^x(2) + Na_in.^x(2)));
% I_NKA_Nahill = @(x, Na_in)  (Na_in./(x(1) + Na_in)).^x(2);
Nahill_fit = lsqcurvefit(I_NKA_Nahill,[10; 1.3;],Na_in_data,I_data);

figure(1);
Na_ins = linspace(0,100);
%normalized
plot(Na_ins,I_NKA_huguet(K_out,Na_ins)./max(abs(I_NKA_huguet(K_out,Na_ins))),...
    'b-','LineWidth',4); hold on
plot(Na_ins,I_NKA_kager(K_out,Na_ins)./max(abs(I_NKA_kager(K_out,Na_ins))),...
    'r--','LineWidth',4);
plot(Na_ins,I_NKA_cressman(K_out,Na_ins)./max(abs(I_NKA_cressman(K_out,Na_ins))),...
    'm:','LineWidth',4);
plot(Na_ins,I_NKA_oschmann(K_out,Na_ins)./max(abs(I_NKA_oschmann(K_out,Na_ins))),...
    '-+','Color',[0.5 0 0.8],'LineWidth',4);
plot(Na_ins,I_NKA_Nahill(Nahill_fit,Na_ins),...
    'k-.','LineWidth',4);
xlabel('[Na^+]_i (mM)'); ylabel('current, normalized');set(gca,'FontSize',20);
legend('Nakao & Gadsby 1989 data', 'Huguet et al. 2016 model', ...
    'Kager et al. 2000 model', 'Cressman et al. 2009 model',...
    'Oschmann et al. 2017 model','proposed NKA model',...
    'Location','southeast');

%data from Figure 11A
[K_out_data, I_data2] = csvimport('nakao_gadsby_1989_fig11a.csv',...
                    'columns', [1, 2] ,'noHeader', true);
%renormalize Idata since its weird
I_data2 = I_data2./max(abs(I_data2));
% Na_out = 150; %mM
Na_in = 50; %mM
% V = 0; %mV

figure(2); clf;
plot(K_out_data, I_data2, 'o', 'MarkerSize',10, 'LineWidth',4); hold on

%my model
I_NKA_Khill = @(x, K_out)  (K_out.^x(2)./(x(1).^x(2) + K_out.^x(2)));
% I_NKA_Khill = @(x, K_out)  (K_out./(x(1) + K_out)).^x(2);
Khill_fit = lsqcurvefit(I_NKA_Khill,[1.5; 1;],K_out_data,I_data2);

                   
figure(2);
K_outs = linspace(0,10);
plot(K_outs,I_NKA_huguet(K_outs,Na_in)./max(abs(I_NKA_huguet(K_outs,Na_in))),...
    'b-','LineWidth',4); hold on;
plot(K_outs,I_NKA_kager(K_outs,Na_in)./max(abs(I_NKA_kager(K_outs,Na_in))),...
    'r--','LineWidth',4);
plot(K_outs,I_NKA_cressman(K_outs,Na_in)./max(abs(I_NKA_cressman(K_outs,Na_in))),...
    'm:','LineWidth',4); 
plot(K_outs,I_NKA_oschmann(K_outs,Na_in)./max(abs(I_NKA_oschmann(K_outs,Na_in))),...
     '-+','Color',[0.5 0 0.8],'LineWidth',4);
plot(K_outs,I_NKA_Khill(Khill_fit,K_outs),...
    'k-.','LineWidth',4);
xlabel('[K^+]_e (mM)'); ylabel('current, normalized');set(gca,'FontSize',20);
legend('Nakao & Gadsby 1989 data', 'Huguet et al. 2016 model', ...
    'Kager et al. 2000 model', 'Cressman et al. 2009 model',...
    'Oschmann et al. 2017 model','proposed NKA model',...
    'Location','southeast');

% %new model only
% I_NKA = @(K_out,Na_in) I_NKA_Nahill(Nahill_fit,Na_in).*I_NKA_Khill(Khill_fit,K_out);
% figure(3);
% plot(Na_in_data, I_data, 'o', 'MarkerSize',10, 'LineWidth',4); hold on
% plot(Na_ins,I_NKA(5.4, Na_ins)./max(abs(I_NKA(5.4, Na_ins))),'k-.','LineWidth',4);
% xlabel('[Na^+]_i (mM)'); ylabel('current, normalized');set(gca,'FontSize',20);
% legend('Nakao & Gadsby 1989 data', 'New model',...
%     'Location','southeast');
% 
% figure(4); 
% plot(K_out_data, I_data2, 'o', 'MarkerSize',10, 'LineWidth',4); hold on
% plot(K_outs,I_NKA(K_outs,50)./max(abs(I_NKA(K_outs,50))),'k-.','LineWidth',4);
% xlabel('[K^+]_e (mM)'); ylabel('current, normalized');set(gca,'FontSize',20);
% legend('Nakao & Gadsby 1989 data', 'New model',...
%     'Location','southeast');


