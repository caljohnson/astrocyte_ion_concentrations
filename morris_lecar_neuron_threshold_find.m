%morris lecar neuron threshold find
%test to see if faster potassium uptake has an effect of neural
%excitability
%finds threshold for the neuron
%plots threshold over time for the resting potential over time traces

addpath('./src'); close all; clear;
load('VKNs.mat');


F = 96485; %C/mol, Faraday's constant
R = 8.31; %J/mol K, ideal gas constant
T = 310; %K, absolute temperature

for kk=1:3
V_Ks = VKNs{kk};
tkk = ts{kk};
I_threshold = size(V_Ks);

for jj=1:size(V_Ks,1)
    
    V_K = V_Ks(jj);

    X0 = [-60.9;0;];
    spiked = 0;
    I_app = 61;
    if kk==3
        I_app = 74;
    end
    tspan = [0,2e2];
    while spiked == 0
        I_app = I_app + 0.2;
        I = @(t) 0.*t + I_app.*(t>10 && t<50);
    
    % options = odeset('RelTol',1e-8,'AbsTol',1e-10,'MaxStep',0.1);
    % [t,X] = ode15s(@(t,X) morris_lecar_neuron_ode(t,X,I,V_K), tspan, X0,options);
        options = odeset('RelTol',1e-6);
        [t,X] = ode45(@(t,X) morris_lecar_neuron_ode(t,X,I,V_K), tspan, X0,options);
    for ii = 1:size(X,1)
        if X(ii,1) > 0
            spiked = 1;
            I_threshold(jj) = I_app;
            break
        end
    end
    end
end

figure(1); hold on
    if kk==1
        plot(V_Ks, I_threshold,'LineWidth',4);
    elseif kk==2
       plot(V_Ks, I_threshold,'--','LineWidth',4);
    else
       plot(V_Ks, I_threshold,'-.','LineWidth',4);
    end

set(gca,'FontSize',20)
xlabel('neuron membrane potential (mV)');
ylabel('threshold current (\muA/cm^2)');

figure(2); hold on
    if kk==1
        plot(tkk,I_threshold,'LineWidth',4);
    elseif kk==2
       plot(tkk, I_threshold,'--','LineWidth',4);
    else
       plot(tkk, I_threshold,'-.','LineWidth',4);
    end

set(gca,'FontSize',20)
xlabel('time (sec)');
ylabel('threshold current (\muA/cm^2)');

end
legend('Elevated potassium, calcium response', ...
    'Elevated potassium, no calcium response',...
    'Resting potassium');
