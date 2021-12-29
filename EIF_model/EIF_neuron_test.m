%EIF neuron test

clear; clc

global V_T del_T g_L g_T E_L C

%paramaters from Fourcaud et al. 2003

C = 1; %muF/cm^2, membrane capacitance
g_L = 0.1;%mS/cm^2, leakage current conductance
g_T = g_L;%threshold current conductance 
del_T = 3.48; %mV, slope factor 
E_L = -65; %mV, resting potential
V_T = -59.9; %mV, soft threshold potential
V_r = -68; %mV, reset potential


%astrocyte-mediated current, reduces excitability
%muA/cm^2
I_ast = @(t) 0.*t;%+ -1.*(t>6e3).*(t<14e3);
% I_ast = @(t) 0.*t + (t-5)./5.*(t>6e3).*(t<14e3);
% I_ast = -0.99; 

%total current - applied + astrocyte-mediated current
%muA/cm^2
I_app = 0.165;
I = @(t) 0.*t + I_app + I_ast(t);

t0 = 0;
tmax = 1e3;
tvec = [t0 tmax];
tmax_final = 5e4;


X0 = V_r;

tic
options = odeset('Events',@depolEvents);
warning('off','all')
for k= 1:200
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

spiked_1 = 0;
spiked_2 = 0;
thresh= 0;
tt = round(size(t,1)/2); %start looking for frequency in second half of simulation trace
while tt <= size(tv,1) && spiked_2 == 0
    if spiked_1 == 1 && xv(tt,1) > thresh && xv(tt-1,1) <=thresh && spiked_2 == 0
        spiked_2 = 1;
        spike_time2 = tt;
    end
    if xv(tt,1) > thresh && xv(tt-1,1) <=thresh && spiked_1 == 0
        spiked_1 = 1;
        spike_time1 = tt;
    end
    tt = tt+1;
end

if spiked_2 == 0
    freq = 0
else
    freq = 1e3.*1/(tv(spike_time2) - tv(spike_time1))
end



function [position,isterminal,direction] = depolEvents(t,y)
global V_T
    if isnan(y)
        position = 0;
    elseif isinf(y)
        position = 0;
    elseif y>=V_T
        position = 0;
    else
        position = 1; % The value that we want to be zero cutoff
    end
    isterminal = 1;  % Halt integration 
    direction = 1;   % The zero can be approached from either direction
end
