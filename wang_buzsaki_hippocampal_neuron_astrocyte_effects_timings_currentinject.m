%wang-buzsaki hippocampal neuron astrocyte effects timingsed current inject
%test to see if faster astrocytic potassium uptake and sodium dynamics have an effect of neural
%excitability

%demonstrates that elevated potassium leads to increased excitability
% and that Ca-transient helps lower excitability to compensate

addpath('./src'); close all; clear;
% load('Na_K_outs.mat');
load('Na_K_outs_superelevatedK.mat');
% load('Na_K_outs_onlyCatransient.mat');

F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature
I = @(t) 0*t + 6.9*(t>5 && t<6); %applied current
% I = @(t) 0.2;
tmax = 1e2;
% tmax = 3e2;

%neural concentrations
K_in = 93.2; %mM  - set to make V_K = -90 mV at rest
Na_in = 17.8; %mM - set to make V_Na = 55 mV at rest

% %neural reversal potentials
% V_Na = (R*T/F).*log(Na_out./Na_in).*1e3;
% V_K = (R*T/F).*log(K_out./K_in).*1e3;


t_checks = [0; 2.5; 5; 7.5; 10; 20; 30;];
% t_checks = [0; 2.5; 5; 7.5; 10; 15; 20; 25; 30; 35;];

%condition 1 = elevated potassium, calcium transient
%condition 2 = elevated potassium, no Ca-transient
%condition 3 = resting potassium, no Ca-transient
for jj=1:3
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
    
    figure(3); plot(V_Nas, V_Ks,'LineWidth',4);  hold on;
    set(gca,'FontSize',20);
    ylabel('neuron potassium reversal potential (mV)');
    xlabel('neuron sodium reversal potential (mV)');

 elseif jj==2
    figure(1); plot(t, V_Nas,'--','LineWidth',4);  hold on;
    ylim([40,60]); set(gca,'FontSize',20);
    xlabel('time (seconds)'); ylabel('neuron sodium reversal potential (mV)');

    figure(2); plot(t, V_Ks,'--','LineWidth',4); hold on;
    ylim([-95,-75]); set(gca,'FontSize',20);
    xlabel('time (seconds)'); ylabel('neuron potassium reversal potential (mV)');
        
    
    figure(3); plot(V_Nas, V_Ks, '--', 'LineWidth',4);  hold on;
    set(gca,'FontSize',20);
    ylabel('neuron potassium reversal potential (mV)');
    xlabel('neuron sodium reversal potential (mV)');
 else
    figure(1); plot(t, V_Nas,'-.','LineWidth',4);  hold on;
    ylim([40,60]); set(gca,'FontSize',20);
    xlabel('time (seconds)'); ylabel('neuron sodium reversal potential (mV)');

    figure(2); plot(t, V_Ks,'-.','LineWidth',4); hold on;
    ylim([-95,-75]); set(gca,'FontSize',20);
    xlabel('time (seconds)'); ylabel('neuron potassium reversal potential (mV)');
    
    figure(3); plot(V_Nas, V_Ks,'-.','LineWidth',4);  hold on;
    set(gca,'FontSize',20);
    ylabel('neuron potassium reversal potential (mV)');
    xlabel('neuron sodium reversal potential (mV)');
    legend('Elevated potassium, calcium response',...
        'Elevated potassium, no calcium response',...
        'Resting potassium/sodium');
 end

    
 
 
end
figure(1); legend('Elevated potassium, calcium response',...
        'Elevated potassium, no calcium response',...
        'Resting potassium/sodium'); xlim([0,40]);
for ii=1:size(t_checks,1)
    plot([t_checks(ii),t_checks(ii)],[-100, 100],'k:','LineWidth',3,'HandleVisibility','off')
end

figure(2); legend('Elevated potassium, calcium response',...
        'Elevated potassium, no calcium response',...
        'Resting potassium/sodium');   xlim([0,40])
for ii=1:size(t_checks,1)
    plot([t_checks(ii),t_checks(ii)],[-100, 100],'k:','LineWidth',3,'HandleVisibility','off')
end

return
V0 = -64;
for ii=1:size(t_checks,1)
for jj=1:3
    
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

    X0 = [V0;0.78;0.088;];
    tvec = [0 tmax];

    [t,X] = ode23(@(t,X) wang_buzsaki_hippocampal_neuron_currentinject_ode(t,X,I,V_Na,V_K), tvec, X0);
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
legend('Elevated potassium, calcium response',...
        'Elevated potassium, no calcium response',...
    'Resting potassium/sodium');
end

