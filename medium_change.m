function [ E_out ] = medium_change(E,n,f,lambda,s_polarization)
%MEDIUM_CHANGE Summary of this function goes here
%   Detailed explanation goes here
F_E=fftshift(fft(E)); % Truncate CO frequencies
sine_i = f/(n/lambda);
sine_t = f/(1/lambda);
rel_inds = abs(sine_t)>1;

cosine_t = sqrt(1-sine_t.^2); % In free space
cosine_i = sqrt(1-sine_i.^2); % Inside crystal

if s_polarization % Fresnel reflection coefficient for each plane wave, from Wiki
    t_s=2*n*cosine_i./(n*cosine_i+cosine_t);
    F_E_out=F_E.*t_s;
else
    t_p=2*n*cosine_i./(cosine_i+n*cosine_t);
    F_E_out=F_E.*t_p;
end
F_E_out(rel_inds) = 0;
E_out = ifft(ifftshift(F_E_out));
end

