%wang-buzsaki hippocampal neuron astrocyte effects timingsed only Ca
%test to see if caclium+sodium dynamics have an effect of neural
%excitability

% demonstrates that Ca-transient lowers excitability

addpath('./src'); close all; clear;
% load('Na_K_outs.mat');
load('Na_K_outs_onlyCatransient.mat');

F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature
I = -0.6; %applied current

%neural concentrations
K_in = 93.2; %mM  - set to make V_K = -90 mV at rest
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest

% %neural reversal potentials
% V_Na = (R*T/F).*log(Na_out./Na_in).*1e3;
% V_K = (R*T/F).*log(K_out./K_in).*1e3;


t_checks = [0; 5;10; 15; 20; 25;30;];

%condition 1 = calcium transient
%condition 2 = no Ca-transient
for jj=1:2
astrocyte_condition = jj;
K_out = K_outs{astrocyte_condition};
Na_out = Na_outs{astrocyte_condition};
t = ts{astrocyte_condition};

%neural reversal potentials
V_Nas = (R*T/F).*log(Na_out./Na_in).*1e3;
V_Ks = (R*T/F).*log(K_out./K_in).*1e3;

 if jj==1
    figure(1); plot(t, V_Nas,'LineWidth',4);  hold on;
    ylim([40,60]); set(gca,'FontSize',20);
    xlabel('time (seconds)'); ylabel('neuron sodium reversal potential (mV)');

    figure(2); plot(t, V_Ks,'LineWidth',4); hold on;
    ylim([-95,-75]); set(gca,'FontSize',20);
    xlabel('time (seconds)'); ylabel('neuron potassium reversal potential (mV)');
 elseif jj==2
    figure(1); plot(t, V_Nas,'--','LineWidth',4);  hold on;
    ylim([40,60]); set(gca,'FontSize',20);
    xlabel('time (seconds)'); ylabel('neuron sodium reversal potential (mV)');

    figure(2); plot(t, V_Ks,'--','LineWidth',4); hold on;
    ylim([-95,-75]); set(gca,'FontSize',20);
    xlabel('time (seconds)'); ylabel('neuron potassium reversal potential (mV)');
 end
   
 
 
end
figure(1); legend('calcium response',...
        'no calcium response'); xlim([0,100]);
for ii=1:size(t_checks,1)
    plot([t_checks(ii),t_checks(ii)],[-100, 100],'k:','LineWidth',3,'HandleVisibility','off')
end
figure(2); legend('calcium response',...
        'no calcium response'); xlim([0,40])
for ii=1:size(t_checks,1)
    plot([t_checks(ii),t_checks(ii)],[-100, 100],'k:','LineWidth',3,'HandleVisibility','off')
end

% V0 = -57.9;
V0 = -55;
% V0s = [-57; -57.5; -58; -60; -65;];
for ii=1:size(t_checks,1)
for jj=1:2
    
    astrocyte_condition = jj;
    K_out = K_outs{astrocyte_condition};
    Na_out = Na_outs{astrocyte_condition};
    t = ts{astrocyte_condition};
    
    ind = find(t==t_checks(ii));
    
    %neural reversal potentials
    V_Nas = (R*T/F).*log(Na_out./Na_in).*1e3;
    V_Ks = (R*T/F).*log(K_out./K_in).*1e3;
    V_Na = V_Nas(ind);
    V_K = V_Ks(ind);

    X0 = [V0;1;0;];
    tvec = [0 1e2];

%     [t,X] = ode23(@(t,X) wang_buzsaki_hippocampal_neuron_komek_version_ode(t,X,I,V_Na,V_K), tvec, X0);
    [t,X] = ode23(@(t,X) wang_buzsaki_hippocampal_neuron_ode(t,X,I,V_Na,V_K), tvec, X0);
    figure(ii+2); hold on
    if jj==1
        plot(t,X(:,1),'LineWidth',4);
    elseif jj==2
        plot(t,X(:,1),'--','LineWidth',4);
    else
        plot(t,X(:,1),'-.','LineWidth',4);
    end 
    set(gca,'FontSize',20);
end
xlabel('time (msec)'); ylabel('neuron membrane potential (mV)');
ylim([-80 40]);
title(['t = ' num2str(t_checks(ii)) ' seconds']);
legend('calcium response',...
        'no calcium response');

end

