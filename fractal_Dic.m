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




%FRACTAL_DIC  computes the degree of intimate contact using the Yang and
%Pitchumani fractal model.
%
%   See also PROCESS_PROFILE_FRACTAL
clear;

%% Surface definition
%according to Y&P
D = 1.32;% 1.2286; % fractal dimension
f = 1.45;% 1.174; % scale factor
h0 = 1.03e-5; % useless, only h0/L0 imports
L0 = h0/0.05;%1.9952e-5;


%with my profile (first shot (1:1278))
D = 1.5162;
h0 = 1.2782e-5;
L0 = 1.0167-4;
f = 1.1504;%1.07539;

%with my profile (everything, filtered)
D = 1.71014;
h0 = 1.31813e-5;
L0 = 2.04736e-4;
f = 1.14277;%1.044;

%cheating ?
%   f = 1.24913;
%   L0=h0/0.03;

s = f ^ ( D/(2-D) );

%% model parameter
Ngen = 15; %number of fractal generation
tfin = 180; %final time
nstep = 180; %number of time step
tspan = linspace(0,tfin,nstep);

%% process parameter
Papp = 276e3; %1340; %pressure applied
T = 350 + 273.16; %temperature

Papp = 1340;
T = 400+273.16;

%% viscosity law
mu = 1.14e-12*exp(26300/T); % viscosity of the fiber matrix bed L&S
%mu = 132.95*exp(2969/T); % viscosity of the fiber matrix bed Mantell and Springer

%% determine useful param for each generation according to Yang and Pitchumani (01)
N = s.^(1:Ngen);
L = (1/f).^(1:Ngen) * L0; % eq (1)
b0 = L./N;

h = (1/f).^((1:Ngen)) * h0; %eq (1)
u = h0 ./ (f.^((1:Ngen+1)-2).*(f+1)); %eq (9)
a0 = [h0 h(1:end-1)] - u(2:end);
af = [h0 h(1:end-1)] - u(1:end-1);

F = L0*Papp./N; %force applied per unit depth on the rectangles

%computing times for squeezing each generation (eq(8))
t = 4*(a0.*b0).^3 / 5 .* ( 1./af.^5  - 1./a0.^5) *mu./F;
% cumulative sum of those times
t(end:-1:1) = cumsum(t(end:-1:1));
t(end+1) = 0;


%% computing Dic
for i=Ngen:-1:1 %for each generation
    indx = (tspan>=t(i+1))&(tspan<t(i)); %which times in tspan are concerned
    if (any(indx)) % if some are,
        Dic(indx) =... compute Dic (eq (12) )
            1/f^i * (5/4*(h0/L0)^2*f^(2*i*D/(2-D)+i+4)/(f+1)^2 *...
            Papp/mu*(tspan(indx)-t(i+1)) + 1).^(1/5);
    end
end

Dic(tspan>t(1)) = 1; %after full squeeze, Dic = 1

%% plotting
plot(tspan, Dic);
hold all;
% add a little circle at each generation change:
plot(t(t<tfin), interp1(tspan,Dic,t(t<tfin)),'o');
