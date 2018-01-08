function out=optodefunUchidaSVU(beta,x)
% function out=optodefunUchidaSVU(beta,x)
% beta - 7 coefficients following Uchida et al. 2008 approach (called 'SVU'
%        in Aanderaa documentation)
% x    - data with temperature as first and phase as second column
%        temperature and phase must have same size
%
% Henry Bittig, GEOMAR/LOV
% 12.05.2017

% out  = (P_0/P_c - 1) / K_SV
% P_0  = c_4 + c_5*T
% P_c  = c_6 + c_7*P
% K_SV = c_1 + c_2*T + c_3*T^2    
% beta(1): c_1
% beta(2): c_2
% beta(3): c_3
% beta(4): c_4
% beta(5): c_5
% beta(6): c_6
% beta(7): c_7


t=x(:,1);
p=x(:,2);

out=((beta(4)+beta(5).*t)./(beta(6)+beta(7).*p) -1)./(beta(1)+beta(2).*t+beta(3).*t.^2);
