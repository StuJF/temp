function hilberts = Hilbert_T(data,dim,n)
% hilberts = Hilbert_T(data,dim,n)
%   Input:
%       data
%       dim: dimension (default=1)
%       n:   number of fft points

if nargin<2
    dim = 1;
end
data = shiftdim(data);
if nargin<3
    n = size(data,dim);
end

fftdata = fft(data,n,dim);
xsize = size(data);
tf = false(size(xsize));
tf(dim) = true;
df = find(~tf);
perm = [find(tf) df];
fftdata = permute(fftdata,perm);
posIdx = 2:floor(n/2)+mod(n,2);
negIdx = ceil(n/2)+1+~mod(n,2):n;
fftdata(posIdx,:) = fftdata(posIdx,:)*2;
fftdata(negIdx,:) = 0;
reperm = [2:find(tf),1,find(tf)+1:length(xsize)];
fftdata = permute(fftdata,reperm);
hilberts = ifft(fftdata,n,dim);