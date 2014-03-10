%PROCESS_PROFILE_FRACTAL  extract the needed Cantor set parameters out of a
%   surface profile.
%
%   [D, h0,L0,f] = PROCESS_PROFILE_FRACTAL(profile) returns the fractal
%   dimension D, the bigest generation rectangle depth h0, the length of
%   the fractal set L0 and the scaling ratio f. profile should be a N-by-2
%   matrix data of a surface profile measurment.  The first column is the
%   abscisse and the second is the depth.
%
%   [D, h0,L0,f] = PROCESS_PROFILE_FRACTAL(profile, plotting) plotting is a
%   binary flag to say whether you want the plots or not the fits.
%
%   [D, h0,L0,f,a0_b0,w0_b0] = PROCESS_PROFILE_FRACTAL(profile) returns,
%   in addition, the geometrical parameters of the first generation
%   recangles (the biggest ones) according to Lee and Springer definition :
%   a0/b0 and w0/b0.
%
%   See also FRACTAL_DIC


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

function [D, h0,L0,f,a0_b0,w0_b0] = process_profile_fractal(profile, plotting)

if (nargin == 1)
    plotting = 1;
end

%% h0: the depth
sig = std(profile(:,2));
h0 = 2*sig;


%% L0 and D: out of the fft
sample_size = profile(end,1)-profile(1,1);
nbplot = size(profile,1);

T = sample_size/nbplot;       % Sample time (resolution)
Fs = 1/T;                     % Sampling frequency
L = nbplot;                   % Length of signal

if (plotting)
    figure(1);
    plot(profile(:,1),profile(:,2));
    title('The profile')
    xlabel('absysse (m)')
end

NFFT = 2^(nextpow2(L)-1); % Next power of 2 from length of profile
Y = fft(profile(:,2),NFFT)/L; % Fourier transform
freq = Fs/2*linspace(0,1,NFFT/2+1);
freq = freq(2:end); %frequency spqn

% Plot single-sided amplitude spectrum.
pow_spec = 2*abs(Y(2:NFFT/2+1)); % power spectrum

if (plotting)
    figure(2)
    loglog(freq,pow_spec,'.');
    title('Power Spectrum of the profile')
    xlabel('Frequency (1/m)')
    ylabel('|Y(1/m)|')
end

%fitting the power law:
Ymodel= @(D,C,L0, Lmin,omega) ...
    C./(omega.^(5-2*D)) .* (omega>1/L0).*(omega<1/Lmin) + ... %power law
    C/((1/L0)^(5-2*D)) * (omega<=1/L0)+... %truncated before L0
    C/((1/Lmin)^(5-2*D)) * (omega>=1/Lmin); % and after Lmin

Yexp = @(omega) interp1(freq,pow_spec,omega);

error = @(D,C,L0,Lmin) norm(Ymodel(D,C,L0,Lmin,freq)-Yexp(freq));

% the actual regression:
options = optimset('MaxFunEvals', 5000, 'MaxIter', 2000);
optim_param = fminsearch(...
    @(param) error(param(1),param(2),param(3),param(4)),...
    [1.2 1e5 3*h0 h0/2],...
    options);

D = optim_param(1);
L0 = optim_param(3);

if (plotting)
    hold all;
    loglog(freq, Ymodel(D,optim_param(2),L0,optim_param(4),freq),'r','LineWidth',3);
end


%% f, the scale factor

Lmodel = @(ff,u) L0*(u+sig)/(ff*h0);
% the famous intercepted length is the number of points above u.
% the intercepted length ratio is then the sum divided by nbplot.
% the intercepeted length per cantor set is that multiplied by L0.
Lexp = @(u) sum(repmat(profile(:,2),1,50) < repmat(u,nbplot,1)) / nbplot * L0;

uspan = linspace(-sig,sig,50);
error = @(ff) norm(Lmodel(ff,uspan)-Lexp(uspan));

if (plotting)
    figure(3);
    plot(uspan,Lexp(uspan),'x');
    hold all;
end
% find f that minimize the error:
f = fminsearch (error, 0.5);

if (plotting)
    plot(uspan, Lmodel(f,uspan));
end

if(nargout>4)
    a0_b0 = h0/L0 * f^((4-D)/(2-D)) / (f+1);
    w0_b0 = (f-1) * f^(D/(2-D));
end

end
