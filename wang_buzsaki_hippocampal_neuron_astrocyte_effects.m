%wang-buzsaki hippocampal neuron astrocyte effects
%test to see if faster astrocytic potassium uptake and sodium dynamics have an effect of neural
%excitability

addpath('./src'); close all; clear;
F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

%neural concentrations
K_in = 93.2; %mM  - set to make V_K = -90 mV at rest
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest

% %neural reversal potentials
% V_Na = (R*T/F).*log(Na_out./Na_in).*1e3;
% V_K = (R*T/F).*log(K_out./K_in).*1e3;

astrocyte_condition = 3;

if astrocyte_condition == 1
    %resting external concentrations
    K_outs = 3.2; %mM 
    Na_outs = 140; %mM
elseif astrocyte_condition == 2
    %astrocyte-effected extracellular concentrations - no calcium transient
    times = [0;     2.73068; 7.27272;  12.812; 30.0656];
    K_outs = [5;        4; 3.41571; 3.245; 3.204;];
    Na_outs = [140; 140.184; 140.294; 140.31; 140.267];
else
    %astrocyte-effected extracellular concentrations - with calcium transient
    times = [0;  3.78391; 5.49047; 7.27272;  8.4968;     9.92;   18.027; 40.5811;   135.913;  324.838;];
    K_outs = [5;    3.82;   3.56;      3.36; 3.28945;  3.20245; 3.15685;  3.18091;  3.20131; 3.19673;];
    Na_outs = [140; 136.117; 132.074; 129.3; 128.722;  129.162; 133.42;  138.04;  140.945; 140.377;];
end


%neural reversal potentials
V_Nas = (R*T/F).*log(Na_outs./Na_in).*1e3;
V_Ks = (R*T/F).*log(K_outs./K_in).*1e3;


% I = 0;% + 0.2.*(t>50 && t<100)+10.*(t>250 && t<300);
I = -0.6;

V0s = [-55; -57; -57.5; -58; -60; -65;];
for jj=1:size(K_outs,1)
for kk=1:size(V0s,1)
    
    V_K = V_Ks(jj);
    V_Na = V_Nas(jj);

    V0 = V0s(kk);
    X0 = [V0;1;0;];
    % X0 = [-10;0;];
    t0 = 0;
    tmax = 1e2;
    tvec = [t0 tmax];

%     [t,X] = ode23(@(t,X) wang_buzsaki_hippocampal_neuron_komek_version_ode(t,X,I,V_Na,V_K), tvec, X0);
    [t,X] = ode23s(@(t,X) wang_buzsaki_hippocampal_neuron_ode(t,X,I,V_Na,V_K), tvec, X0);
    figure(jj); hold on
    title(['[Na^+]_e = ' num2str(Na_outs(jj)) ' mM, [K^+]_e = ' num2str(K_outs(jj)) ' mM']);
    plot(t,X(:,1),'LineWidth',4); ylim([-100,40])   
    set(gca,'FontSize',20);
xlabel('time (seconds)'); ylabel('neuron membrane potential (mV)');
lbls{kk} = ['V_0 = ' num2str(V0) ' mV'];
end
legend(lbls)
% legend('Elevated potassium, no calcium response',...
%         'Elevated potassium, calcium response',...
%         'Elevated potassium, no astrocyte response',...
%     'Resting potassium/sodium');
% legend('Default E_{Na},E_{K}','Lower (restored) E_{Na},E_{K}','Lowest (pathological) E_{Na},E_{K}')


end
% legend(lbls)
% %pulse 1
% plot([50,100],[-90,-90],':k','LineWidth',4);
% plot([50, 50], [-100, -90],':k','LineWidth',4)
% plot([100, 100], [-100, -90],':k','LineWidth',4)
% 
% %pulse 2
% plot([250,300],[-90,-90],':k','LineWidth',4);
% plot([250, 250], [-100, -90],':k','LineWidth',4)
% plot([300, 300], [-100, -90],':k','LineWidth',4)
