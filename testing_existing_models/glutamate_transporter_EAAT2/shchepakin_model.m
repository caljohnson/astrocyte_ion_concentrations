%EAAT model - Shchepakin et al. 2019 model


%smallest system 4-state model
%glutamate present
%EAAT2 coefficients - medians from table 1
m_1p = 1155.84;
m_1n = 101.99;
m_2p = 556.67;
m_2n = 697.34;
m_3p = 422.36;
m_4p = 31.12;
m_4n = 13.35;
%current coefficients - medians from table 1
A = 167.29;
B = 265.48;
D = 163.93;

%coefficient matrix w/ glutamate present
coeff_matrix_gp = [-m_1p - m_4n, m_1n, 0, m_4p; ...
                m_1p, -m_1n - m_2p, m_2n, 0; ...
                0, m_2p, -m_2n - m_3p, 0; ...
                m_4n, 0, m_3p, -m_4p;];
            
%w/ no glutamate present
m_1p = 0;
coeff_matrix_NOgp = [-m_1p - m_4n, m_1n, 0, m_4p; ...
                m_1p, -m_1n - m_2p, m_2n, 0; ...
                0, m_2p, -m_2n - m_3p, 0; ...
                m_4n, 0, m_3p, -m_4p;];
                                    
%ODE solving stuff            
ode_rhs_gp = @(t,X) coeff_matrix_gp*X;
ode_rhs_NOgp = @(t,X) coeff_matrix_NOgp*X;

x0 = [m_4p/(m_4p+m_4n); 0; 0; m_4n/(m_4p+m_4n);];

%glutamate pulse parameters
T_pulse1 = 0.05; %first Glu pulse, 50 ms
T_pulse2 = 0.03; %second Glu pulse, 30 ms
%variable delay between pulses
T_delay = [0.005; 0.010; 0.015; 0.020; 0.030; 0.040; ...
           0.050; 0.060; 0.080; 0.100; 0.150; 0.200; 0.250; 0.300;];
TF = 0.5;

for ii = 1:size(T_delay,1)

[t1,x1] = ode23s(ode_rhs_gp,[0,T_pulse1],x0);
[t2,x2] = ode23s(ode_rhs_NOgp,[t1(end),t1(end)+T_delay(ii)],x1(end,:));
[t3,x3] = ode23s(ode_rhs_gp,[t2(end),t2(end)+T_pulse2],x2(end,:));
[t4,x4] = ode23s(ode_rhs_NOgp,[t3(end),TF],x3(end,:));

%convert system states to transporter current
t = [t1;t2;t3;t4;];
x = [x1;x2;x3;x4;];

I = -A*x(:,1) - B*x(:,2) - D*x(:,4) + (A*m_4p +D*m_4n)/(m_4p + m_4n);
plot(t,I,'LineWidth',2); hold on
end


            
            