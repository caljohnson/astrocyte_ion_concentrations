%wang-buzsaki hippocampal neuron test
%test to see if faster potassium uptake has an effect of neural
%excitability

addpath('./src'); close all; clear;
% load('VKNs.mat');


F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

% V_Ks = [-83.5; -77.5; -85.5;]
% V_Ks = [-84.5;-87; -78; -90;];
% V_Nas = [54.1; 46.35; 55; 55;];
% V_Ks = -90; %from paper
% V_Ks = [-90;-78;] %paper, elevated potassium!
% V_Nas = [55; 55;] %from paper
% V_Ks = -55;
V_Ks = -90; %resting
V_Nas = 55; %resting
% V_Ks = -75;
% V_Nas = 53.5;
V_Nas = 55;
tmax = 1e3;

% I = 5;
I = 0.11;
% I = 0;

V0s = -50;
% V0s = [-58.8; -58.65; -58.5; -58;];
for kk=1:size(V0s,1)
for jj=1:size(V_Ks,1)
    
    V_K = V_Ks(jj);
    V_Na = V_Nas(jj);

    V0 = V0s(kk);
%     X0 = [V0;1;0;];
%     X0 = [-55;1;0;];
    X0 = [V0;0.78;0.088;];
    
    
    tvec = [0 tmax];

%     [t,X] = ode23s(@(t,X) wang_buzsaki_hippocampal_neuron_komek_version_ode(t,X,I,V_Na,V_K), tvec, X0);
    [t,X] = ode23s(@(t,X) wang_buzsaki_hippocampal_neuron_ode(t,X,I,V_Na,V_K), tvec, X0);
    figure(jj); hold on
%     figure(1); hold on
    plot(t,X(:,1),'LineWidth',4); ylim([-100,40])   
    set(gca,'FontSize',20);
xlabel('time (msec)'); ylabel('neuron membrane potential (mV)');

    figure(jj+2);
    plot(X(:,3), X(:,2),'o','LineWidth',4);
    set(gca,'FontSize',20);
    xlabel('n'); ylabel('h');
    
    figure(jj+4); hold on
    plot(t,X(:,2),'LineWidth',4);
    plot(t,X(:,3),'LineWidth',4);
    set(gca,'FontSize',20);
    xlabel('time (msec)'); ylabel('gating variable');
    legend('h','n');
end
% legend('Elevated potassium, no calcium response',...
%         'Elevated potassium, calcium response',...
%         'Elevated potassium, no astrocyte response',...
%     'Resting potassium/sodium');
% legend('Default E_{Na},E_{K}','Lower (restored) E_{Na},E_{K}','Lowest (pathological) E_{Na},E_{K}')
% lbls{kk} = ['V_0 = ' num2str(V0) ' mV'];

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

% 
%see if neuron can "turn off" after astrocyte effects
% V_Na = 54;
% V_K = -70;
I= -0.8;
[t2,X2] = ode23s(@(t,X) wang_buzsaki_hippocampal_neuron_ode(t,X,I,V_Na,V_K), [tmax, 2*tmax], X(end,:));
t = [t; t2;];
X = [X; X2;];
figure(6);
plot(t,X(:,1),'LineWidth',4); ylim([-80,40]);
% xlim([680 730]);
set(gca,'FontSize',20);
xlabel('time (msec)'); ylabel('neuron membrane potential (mV)');


spiked_1 = 0;
spiked_2 = 0;
tt = 1;
thresh = -20;
while tt <= size(t,1) && spiked_2 == 0
    if spiked_1 == 1 && X(tt,1) > thresh && X(tt-1) <=thresh && spiked_2 == 0
        spiked_2 = 1;
        spike_time2 = tt;
    end
    if X(tt,1) > thresh && X(tt-1) <=thresh && spiked_1 == 0
        spiked_1 = 1;
        spike_time1 = tt;
    end
    tt = tt+1;
end
    
if spiked_1 == 0
    delay = NaN;
else
    delay = t(spike_time1);
end

if spiked_2 == 0
    freq = 0;
else
    freq = 1e3.*1/(t(spike_time2) - t(spike_time1));
end

freq

