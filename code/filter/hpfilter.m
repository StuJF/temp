function Hd = hpfilter(Fc,Fs,N)
% butterworth highpass filter
if nargin<3
    N=10;% Order
end
  
Fc = 0.5;  % Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.highpass('N,F3dB', N, Fc, Fs);
Hd = design(h, 'butter');

% [EOF]
