%morris lecar neuron test
%test to see if faster potassium uptake has an effect of neural
%excitability

addpath('./src'); close all; clear;
load('VKNs.mat');


F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

V_Ks = [-83.5; -77.5; -85.5;]

% I = @(t) 0.*t + 73.*(t>10);
% I = @(t) 0.*t + 73.*(t>10 && t<50);
I = @(t) 0.*t + 80.*(t>10 && t<50);
% I = @(t) 0.*t + 0.*(t>10 && t<50);

for jj=1:3
    
    V_K = V_Ks(jj);

    X0 = [-60.9;0;];
    % X0 = [-10;0;];
    tspan = [0,2e2];
    % options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.1);
    % [t,X] = ode15s(@(t,X) morris_lecar_neuron_ode(t,X,I,V_K), tspan, X0,options);
    options = odeset('RelTol',1e-6);
    [t,X] = ode45(@(t,X) morris_lecar_neuron_ode(t,X,I,V_K), tspan, X0,options);

    figure(1); hold on
    if jj==1
        plot(t,X(:,1),'LineWidth',4); ylim([-100,40])
    elseif jj==2
        plot(t,X(:,1),'--','LineWidth',4); ylim([-100,40])
    else
        plot(t,X(:,1),'-.','LineWidth',4); ylim([-100,40])
    end
end

plot([10,50],[-90,-90],':k','LineWidth',4);
plot([10, 10], [-100, -90],':k','LineWidth',4)
plot([50, 50], [-100, -90],':k','LineWidth',4)
legend('Elevated potassium, calcium response', ...
    'Elevated potassium, no calcium response',...
    'Resting potassium',...
    'Applied current timing, 73 \muA/cm^2');
set(gca,'FontSize',20)
xlabel('time (msec)'); ylabel('neuron membrane potential (mV)');

return
V_Ks = [-88.7; -81.9; -85.5;];

% I = @(t) 0.*t + 75.*(t>10);
I = @(t) 0.*t + 76.*(t>10 && t<50);
for jj=1:3
    
    V_K = V_Ks(jj);

    X0 = [-60.9;0;];
    % X0 = [-10;0;];
    tspan = [0,2e2];
    % options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.1);
    % [t,X] = ode15s(@(t,X) morris_lecar_neuron_ode(t,X,I,V_K), tspan, X0,options);
    options = odeset('RelTol',1e-6);
    [t,X] = ode45(@(t,X) morris_lecar_neuron_ode(t,X,I,V_K), tspan, X0,options);

    figure(2); hold on
    plot(t,X(:,1),'LineWidth',4); ylim([-100,40])
end
legend('Elevated potassium, calcium response', ...
    'Elevated potassium, no calcium response',...
    'Resting potassium');
set(gca,'FontSize',20);
xlabel('time (seconds)'); ylabel('neuron membrane potential (mV)');

