%INTIMATE_CONTACT  computes the degree of intimate contact using the Lee
%and Springer model (isothermal)
%
%   Dic = INTIMATE_CONTACT(P,Tc,t) returns the degree of intimate contact
%   Dic at time t (s) for a given pressure P(Pa) and a given temperature
%   T(c). note that t can either be a scalar or a vector of times.The Lee
%   and Springer values for the viscosity Arrhenius law is used.

%This file is part of Fractal Dic.
%
%    Fractal Dic is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    Fractal Dic is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with ATP Simulation.  If not, see <http://www.gnu.org/licenses/>.
%
%
% Copyright 2011 Arthur Levy


function Dic = intimate_contact (P,TC, t)

ratio = 1;
w0_b0 = ratio* 0.9256 + (ratio-1)*ratio;%0.74;
a0_b0 = ratio^2*0.0883;%1.23;

w0_b0 = 1;%24;%1;%1;%0.743;
a0_b0 = 9.69e-3;%1.7;%0.3;%1.237;

T = TC+273.16;
%Lee & Springer 87
%  mumf = 1.14 * 10^-12 *...
%     exp(26300/T);

% mumf = 1.13 * 10^-10 *...
%     exp(19100/T);
mumf = 132.95*exp(2969/T);
% mumf = 643*exp(4367/T);


Dic = 1 / ( 1 + w0_b0) * ...
    ( 1 + ...
    5*P / mumf * (1 + w0_b0) * (a0_b0)^2 * t ) .^ (1/5);

Dic(Dic>1)=1;
