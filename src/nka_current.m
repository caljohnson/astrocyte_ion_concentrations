function [I_current] = nka_current(Na_in,K_out,type)
%NKA_CURRENT model of the NKA (Sodium-Potassium Pump) current
%   relies on concentrations of internal sodium (Na_in) and external potassium (K_out) 
%   Uses either the NKA model by Kager et al. 2000 (type = 1)
%                                Huguet et al. 2016 (type = 2)
%                             or Cressman et al. 2009 (type = 3)

if nargin < 4 || isempty(type)
    type = 4;
end

if(type==1)%Kager et al. 2000
    kK = 3.5; %mM, half activation [K+]e concentration for the astrocytic NKA
    kNa = 10; %mM, half activation [Na+]i concentration for the astrocytic NKA
    I_current = (K_out./(kK + K_out)).^2.*(Na_in./(kNa + Na_in)).^3;
    
elseif(type==2) %Huguet et al. 2016
    kK = 2; %mM, half activation [K+]e concentration for the astrocytic NKA
    kNa = 7.7; %mM, half activation [Na+]i concentration for the astrocytic NKA
    I_current = (K_out./(kK + K_out)).^2.*(Na_in./(kNa + Na_in)).^3;
    
elseif(type==3) %Cressman et al. 2009
    I_current = (1+exp((25-Na_in)/3)).^(-1).*(1+exp(5.5-K_out)).^(-1);% from paper
    
elseif(type==4) %my own, based on Nakao & Gadsby 1989 data
    I_NKA_Nahill = @(x, Na_in)  (Na_in.^x(2)./(x(1).^x(2) + Na_in.^x(2)));
    I_NKA_Khill = @(x, K_out)  (K_out.^x(2)./(x(1).^x(2) + K_out.^x(2)));
    Nahill_fit = [10.5; 1.3;];
    Khill_fit = [1.1; 1.15;];
    I_current = I_NKA_Nahill(Nahill_fit,Na_in).*I_NKA_Khill(Khill_fit,K_out);
end

% K_flux = -2*I_current;
% Na_flux = 3*I_current;

end

