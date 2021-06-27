function [signal,zi] = fir_filt(b,signal,zi)
% filter data using fir filter
% [signal,zi] = fir_filt(b,signal,zi)
%  input:
%       - b:    fir filter coefficent
%       -signal:
%       -zi:    initial condition
if nargin<3
    zi=[];
end


if mod(b,2)==0
    error('filter order not even')
end


blocksize=1000;%max number of chunk for filt

[~,nsam]=size(signal);

delay=length(b);

if isempty(zi)
    signal = [repmat(2*signal(:,1),1,delay) - signal(:,1+mod(((delay+1):-1:2)-1,size(signal,2))),signal];
    
    cut=true;
    newsam=nsam+length(b);
else
    cut=false;
    newsam=nsam;
end



fnum=ceil(newsam/blocksize);%split signal to chunk
for id=1:fnum
    
    range=(1+(id-1)*blocksize):min(id*blocksize,newsam);
    [signal(:,range),zi]=filter(b,1,signal(:,range),zi,2);
    
end


if cut
    signal=signal(:,end-nsam+1:end);
end


end

