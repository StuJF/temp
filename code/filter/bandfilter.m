function hd=bandfilter(w,fs,N)
% butterworth bandpass filter
% hd=bandfilter(w,fs,N)
if nargin<3
    N=10;% Order
end
Fs = fs;  % Sampling Frequency

 
Fc1 = w(1);  % First Cutoff Frequency
Fc2 = w(2);  % Second Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
hd = design(h, 'butter');
end

