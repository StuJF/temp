function b = fir_filter(w,fs)
%fir_filter Returns a discrete-time fir filter b using firls
%fir_filter(w,fs)
% - w is the frequency in HZ;[w1 w2] 
% - fs is sample frequency in HZ
% - lowpass:[0 high] highpass:[low 0] bandpass:[low high]
if nargin<2
    error('fir_filter require sampling frequency');
end

lowcutoff = w(1);
highcutoff = w(2);
min_filter_order = 15;
trans = 0.15;
min_fac = 3;
if lowcutoff>0
    filtorder = min_fac*ceil(fs/lowcutoff);
elseif highcutoff>0
    filtorder = min_fac*ceil(fs/highcutoff);
end
if filtorder < min_filter_order
    filtorder = min_filter_order;
end
nqy = fs/2;
if highcutoff>nqy
    error('high cutoff>fs/2');
end
if lowcutoff>nqy
    error('low cutoff>fs/2');
end
if highcutoff*(1+trans)>nqy
   error('high cutoff frequency too close to Nyquist frequency'); 
end
if lowcutoff<0
    error('low cutoff<0');
end
if lowcutoff>0
    if highcutoff>0 % bandpass
        f = [0 lowcutoff*(1-trans)/nqy lowcutoff/nqy highcutoff/nqy highcutoff*(1+trans)/nqy 1];
        m = [0 0 1 1 0 0];
    else % highpass
        f = [0 lowcutoff*(1-trans)/nqy lowcutoff/nqy 1];
        m = [0 0 1 1];
    end
elseif highcutoff>0 % lowpass
    f = [0 highcutoff/nqy highcutoff*(1+trans)/nqy 1];
    m = [1 1 0 0];
else
    error('at least provide one non-zero low cutoff or high cutoff')
end
    
b = firls(filtorder,f,m);







end



% [EOF]
