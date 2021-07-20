%
%  Morris-Lecar model of an excitable cell (barnacle muscle)
%             Victor Matveev, June 10 2010
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
figure('Position',[1 200 900 350]);  % position the window

global C V1 V2 V3 V4 gCa gK gL VK VL VCa Iapp vv;

V1  = -1.2; V2 = 18;  V3 = 2; V4 = 30;
gCa = 4.4;  gK = 8;   gL = 2;
VCa = 120;  VK = -84; VL=-60;
C = 20; Iapp = 0;

minf = @(V) (1 + tanh((V-V1)/V2) )/2;
winf = @(V) (1 + tanh((V-V3)/V4) )/2;
r    = @(V) cosh((V-V3)/(2*V4)) / 25;

Itotal = @(V,w) -gCa*minf(V).*(V-VCa)-gK*w.*(V-VK)-gL*(V-VL)+Iapp;
f = @(t,y) [ Itotal(y(1),y(2))/C; r(y(1)).*(winf(y(1))-y(2)) ];
          
options = odeset('RelTol',1e-6);

T = 200;
IC = [-10 0];
[T Y] = ode45(f, [0 T], IC);

subplot(1,2,1);          % Time plot
plot(T, Y(:,1), 'r-');
xlabel('time'); ylabel('V(t)');

subplot(1,2,2);          % Phase plane plot
plot(Y(:,1), Y(:,2), 'r-');
hold on;

i = 1;
vRange = -75:5:75;

for vv = vRange   % This calculates the V-nullcline
 w(i) = fsolve(@(w) Itotal(vv,w), 0);
 i = i + 1;
end;

plot(vRange, w,'k-');              % plot the V-nullcline
plot(vRange, winf(vRange),'b-');   % Plot the w-nullcline
xlabel('V(t)'); ylabel('w(t)'); title('Phase space');
axis([-75 75 -0.1 0.5]);
legend('Trajectory', 'V-nullcline', 'w-nullcline', 'Location', 'SE');