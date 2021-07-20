%izhikevich hippocampal neuron test
%test to see if faster potassium uptake has an effect of neural
%excitability

addpath('./src'); close all; clear;
% load('VKNs.mat');


F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

% V_Ks = [-83.5; -77.5; -85.5;]
V_Ks = [-90;-90; -90;];

I = @(t) 0.*t + 100.*(t>50 && t<100)+100.*(t>250 && t<300);
for jj=1:3
    
    V_K = V_Ks(jj);

    X0 = [-60.9;0;];
    % X0 = [-10;0;];
    t0 = 0;
    tmax = 1e3;
    tvec = [t0 tmax];

    options = odeset('Events',@depolEvents);
    for k= 1:7
        [t{k},X{k}] = ode45(@(t,X) izhikevich_hippocampal_neuron_ode(t,X,I,V_K), tvec, X0,options);
        X0 = [-50; 50;];
        tvec = [t{k}(end) t{k}(end)+tmax];
    end
    tv = cat(1,t{:});
    xv = cat(1,X{:});

    figure(2); hold on
    plot(tv,xv(:,1),'LineWidth',4); ylim([-100,40])
end
% legend('Elevated potassium, calcium response', ...
%     'Elevated potassium, no calcium response',...
%     'Resting potassium');
set(gca,'FontSize',20);
xlabel('time (seconds)'); ylabel('neuron membrane potential (mV)');

%pulse 1
plot([50,100],[-90,-90],':k','LineWidth',4);
plot([50, 50], [-100, -90],':k','LineWidth',4)
plot([100, 100], [-100, -90],':k','LineWidth',4)

%pulse 2
plot([250,300],[-90,-90],':k','LineWidth',4);
plot([250, 250], [-100, -90],':k','LineWidth',4)
plot([300, 300], [-100, -90],':k','LineWidth',4)


function [position,isterminal,direction] = depolEvents(t,y)
    position = y(1) - 40; % The value that we want to be zero 40mV cutoff
    isterminal = 1;  % Halt integration 
    direction = 1;   % The zero can be approached from either direction
end
