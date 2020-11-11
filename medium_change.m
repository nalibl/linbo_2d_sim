function [ co_E_interp ] = medium_change(E,n,f,s_polarization)
%MEDIUM_CHANGE Summary of this function goes here
%   Detailed explanation goes here
co_freq=max(abs(f))/n; % Find cut-off frequency for birefringent medium
co_lim=abs(f)<co_freq; % Get the indices of the angles below the CO
F_E=fftshift(fft(E)); % FFT of signal at end of crystal, inside
F_E=F_E(co_lim); % Truncate CO frequencies
co_f=f(co_lim)/max(abs(f(co_lim)))*max(abs(f)); % Renormalize to new frequencies
co_f_n=co_f/max(abs(co_f)); % Sine of angle, normalized frequency

cosine_t = sqrt(1-co_f_n.^2); % In free space
cosine_i = sqrt(1-(co_f_n/n).^2); % Inside crystal

if s_polarization % Fresnel reflection coefficient for each plane wave
    t_s=2*n*cosine_i./(n*cosine_i+cosine_t);
    F_E=F_E.*t_s;
else
    t_p=2*n*cosine_i./(cosine_i+n*cosine_t);
    F_E=F_E.*t_p;
end
co_E = ifft(ifftshift(F_E)); % Spatial domain after medium change
df = f(2)-f(1);
L = 1/df;
dx_old = 1/(2*max(abs(f)));
dx_new = 1/(2*max(abs(f(co_lim))));%1/(co_f(end)-co_f(1));
y_new = -L/2:dx_new:L/2;
y_old = -L/2:dx_old:L/2-dx_old;
co_E_interp = interp1(y_new, co_E, y_old);
end

