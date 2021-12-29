%EIF neuron test

clear; clc

global V_T del_T g_L g_T E_L

g_L = 0.1;%leakage current conductance
g_T = g_L;%threshold current conductance - adjust to account for astrocyte effects
del_T = 3.6; %slope factor - adjust to account for astrocyte effects
E_L = -70; %mV, resting potential

%astrocyte-mediated current, reduces excitability
I_ast = @(t) 0.*t + -1.*(t>6e3).*(t<14e3);
% I_ast = @(t) 0.*t + (t-5)./5.*(t>6e3).*(t<14e3);
% I_ast = -0.99; 

%total current - applied + astrocyte-mediated current
I = @(t) 0.*t + 3 + I_ast(t);


t0 = 0;
tmax = 1e3;
tvec = [t0 tmax];
tmax_final = 5e4;

V_T = -46; %soft threshold potential
V_r = -60; %reset potential

X0 = V_r;

tic
options = odeset('Events',@depolEvents);
warning('off','all')
for k= 1:2000
    [t{k},X{k}] = ode23s(@(t,X) EIF_neuron_ode(t,X,I), tvec, X0,options);
    if t{k}(end) == tvec(end)
        X0 = X{k}(end);
    else
        X0 = V_r;
    end
    if t{k}(end) > tmax_final
        break
    else
        tvec = [t{k}(end) t{k}(end)+tmax];
    end
end
toc
tv = cat(1,t{:});
xv = cat(1,X{:});

figure(2); clf; 
plot(tv.*1e-3,xv(:,1),'LineWidth',4); ylim([-100,40])
set(gca,'FontSize',20);
xlabel('time (seconds)'); ylabel('neuron membrane potential (mV)');

% %pulse 1
% plot([50,100],[-90,-90],':k','LineWidth',4);
% plot([50, 50], [-100, -90],':k','LineWidth',4)
% plot([100, 100], [-100, -90],':k','LineWidth',4)
% 
% %pulse 2
% plot([250,300],[-90,-90],':k','LineWidth',4);
% plot([250, 250], [-100, -90],':k','LineWidth',4)
% plot([300, 300], [-100, -90],':k','LineWidth',4)


function [position,isterminal,direction] = depolEvents(t,y)
    if isnan(y)
        position = 0;
    elseif isinf(y)
        position = 0;
    else
        position = 1; % The value that we want to be zero cutoff
    end
    isterminal = 1;  % Halt integration 
    direction = 1;   % The zero can be approached from either direction
end
