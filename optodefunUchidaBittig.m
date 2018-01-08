function out=optodefunUchidaBittig(beta,x)
% function out=optodefunUchidaBittig(beta,x)
% beta - 7 coefficients following Bittig et al. 2017 approach (modification
%        of Uchida et al. 2010 equation) 
% x    - data with temperature as first and phase as second column
%        temperature and phase must have same size
%
% Henry Bittig, GEOMAR/LOV
% 12.05.2017

% out  = ((1 + c_4*T)/(c_5 + c_6*P + c_7*P^2) - 1) / (c_1 + c_2*T + c_3*T^2)
% beta(1): c_1
% beta(2): c_2
% beta(3): c_3
% beta(4): c_4
% beta(5): c_5
% beta(6): c_6
% beta(7): c_7


t=x(:,1);
p=x(:,2);

out=((1+beta(4).*t)./(beta(5)+beta(6).*p+beta(7).*p.^2) -1)./(beta(1)+beta(2).*t+beta(3).*t.^2);
