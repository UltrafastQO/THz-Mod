function [Ef,f,r,phi,out_slope] = THz_fft(Et,t,n_pad,f_fit,f_r,in_slope)

% This function performs the FFT of a THz field transient acquired via
% electro-optic sampling
% Author: Nicolas Couture, 2019

% Et       : EOS Signal
% t        : EOS delay in ps
% n_pad    : number of zeros to add to the end of the transient
% f_fit    : The phase should be extrapolated from f_fit to zero (optional)
% f_r      :  Right limit of the data to be used for extrapolation. (optional)
% fitslope : given slope with which the phase should be fitted. (optional)
%
% returns
% Ef        : Complex FFT
% f         : Frequency (THz)
% r         : amplitude
% phi       : phase
% out_slope : Fit_slope

global Ef_0;

% Check input parameters
if nargin < 3
    n_pad =0;
end

if nargin < 4
  f_fit = 0;
  f_r = 0;
end

% set flag, when the slope of the phase shall be fit by external given slope.
if (nargin==6) slope_ext = 1; else slope_ext = 0;  out_slope = 0; end



%sampling frequency in [THz]
timestep = mean(diff(t));           % Timestep in [ps]
Fs = 1./timestep;                   % Sampling frequency in [THz]

% Length of Ef is calculated from the length of Et 
% and the desired number of appended zeros, so n_FFT is a square number
% 2^(nextpow2(length(Et)+n_pad)); commented out
n_FFT = length(Et)+n_pad;         % Number of amplitudes to be calculated.

%Subtrahiere Offset
%Et_womean= Et-mean(Et);

Ef_0 = fft(Et,n_FFT);                % Fast Fourier Transformation

NumUniquePts = ceil((n_FFT+1)/2);   % Number of non-redundant data points
                                    % Output from FFT has the form:
                                    % [A  B C D  E E  D C B  ] if n_FFT odd
                                    % [A  B C D   E   D C B  ] if n_FFT even
%Reduce Fourier spectra to non-redundant data points because 
% amplitudes for negative and positive frequencies are identical.
Ef = Ef_0(1:NumUniquePts);

% Pre-assignment of various variables
f = zeros(NumUniquePts,1);          % frequency-matrix
phi = zeros(NumUniquePts,1);        % Phase
delta_phi = zeros(NumUniquePts,1);  % delta_phi @@@ 
%Frequency scale
f = (0:NumUniquePts-1)*Fs/n_FFT;    % Occupancy of the frequency vector.
f = f';

% Normalization of the amplitude according to “Plancherel’s theorem”: int |f(w)|^2 dw = int |f(t)|^2 dt
r = abs(Ef);

intt = timestep*sum(Et.^2);
intf = mean(diff(f))*sum(r.^2);
if (intf~=0) a1 = sqrt(intt/intf); else a1 = 1; end
r = r.*a1;

if ~rem(n_FFT,2)                    % ???
    r(end) = r(end)/2;              % ???
end                                 % ???


% Calculation of Phase
phi = unwrap(angle(Ef));%-unwrap(angle(Ef));  %Minus sign gives corrected engineering convention 

% Correct phase for low frequencies by extrapolation.
if (f_fit > f(1))
  % Shift coordinates so that f=f_fit is at the origin, with phi(f_fit) = 0;
  f_fitrange = f(find(( f>f_fit) & (f<f_r)));                              % Frequency range to be extrapolated.
  f_offset = f_fitrange(1);                                                % So that extrapolation can be carried out continuously
  f_fitrange = f_fitrange - f_offset;                                      % - " -

  phi_fitrange = phi(find(( f>f_fit) & (f<f_r)));                          % Phase in the corresponding frequency range.
  phi_offset = phi_fitrange(1);                                            % So that extrapolation can be carried out continuously
  phi_fitrange = phi_fitrange - phi_offset;                                % - " -
%  phi_fit = fit(f_fitrange,phi_fitrange,'a.*x + b.*x.^2','StartPoint',[0 0]); % search fit parameter
  if (slope_ext == 0)
    phi_fit  = fit(f_fitrange,phi_fitrange,'a.*x','StartPoint',[0]);          % search fit parameter
    out_slope = phi_fit;
  else phi_fit = in_slope;
  end
  f_extra = f(1:find(f>f_fit,1,'first'))-f_offset;                         % generate frequencies for extrapolation
%  phi(1:find(f>f_fit,1,'first')) = phi_offset + phi_fit.a.*f_extra -
%  0.01.*phi_fit.b.*f_extra.^2; % extrapolate
  phi(1:find(f>f_fit,1,'first')) = phi_offset + phi_fit.a.*f_extra;        % extrapolate
  phi = phi - phi(1);                                                      % correct offset
end


Ef = r.*exp(1i.*phi);

