%morris lecar neuron test
%test to see if faster potassium uptake has an effect of neural
%excitability

addpath('./src'); close all; clear;
load('VKNs.mat');


F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

V_Ks = [-83.5; -77.5; -85.5;]

%---- Three pulses - vary strength

% I = @(t) 0.*t + 73.*(t>10);
% I_app = 180;
I_app = 120;
I = @(t) 0.*t + I_app.*(t>10 && t<20) + I_app.*(t>35 && t<45) + I_app.*(t>60 && t<70);
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

%pulse 1
plot([10,20],[-90,-90],':k','LineWidth',4);
plot([10, 10], [-100, -90],':k','LineWidth',4)
plot([20, 20], [-100, -90],':k','LineWidth',4)

%pulse 2
plot([35,45],[-90,-90],':k','LineWidth',4);
plot([35, 35], [-100, -90],':k','LineWidth',4)
plot([45, 45], [-100, -90],':k','LineWidth',4)

%pulse 3
plot([60,70],[-90,-90],':k','LineWidth',4);
plot([60, 60], [-100, -90],':k','LineWidth',4)
plot([70, 70], [-100, -90],':k','LineWidth',4)

legend('Elevated potassium, calcium response', ...
    'Elevated potassium, no calcium response',...
    'Resting potassium',...
    ['Applied current timing, ' num2str(I_app) ' \muA/cm^2']);
set(gca,'FontSize',20)
xlabel('time (msec)'); ylabel('neuron membrane potential (mV)');


%---- Two pulses - vary timing
% V_Ks = [-83.5; -77.5; -85.5;]

I_app = 180;
t2 = 55;
I = @(t) 0.*t + I_app.*(t>10 && t<20) + I_app.*(t>t2 && t<t2+10);

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
    if jj==1
        plot(t,X(:,1),'LineWidth',4); ylim([-100,40])
    elseif jj==2
        plot(t,X(:,1),'--','LineWidth',4); ylim([-100,40])
    else
        plot(t,X(:,1),'-.','LineWidth',4); ylim([-100,40])
    end
end

%pulse 1
plot([10,20],[-90,-90],':k','LineWidth',4);
plot([10, 10], [-100, -90],':k','LineWidth',4)
plot([20, 20], [-100, -90],':k','LineWidth',4)

%pulse 2
plot([t2,t2+10],[-90,-90],':k','LineWidth',4);
plot([t2, t2], [-100, -90],':k','LineWidth',4)
plot([t2+10, t2+10], [-100, -90],':k','LineWidth',4)


legend('Elevated potassium, calcium response', ...
    'Elevated potassium, no calcium response',...
    'Resting potassium',...
    ['Applied current timing, ' num2str(I_app) ' \muA/cm^2']);
set(gca,'FontSize',20)
xlabel('time (msec)'); ylabel('neuron membrane potential (mV)');
