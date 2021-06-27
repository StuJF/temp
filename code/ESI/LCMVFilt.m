function Wimaging = LCMVFilt(L,Covy,nsource,type,varargin)
% W = LCMVFilt(L,lamda,Covy,type,nsource)
% Minimum-norm method for EEG Source Imaging
% Input:
%   L: Gain matrix
%   Covy: data covariance (default is L*L')
%   lamda: regularization parameters (default is 1/SNR^2)
%   type:  method (default is unit-gain)

[nchan,~] = size(L);
if nargin<2||isempty(Covy)
    Covy = L*L';
end
if nargin<3||isempty(nsource)
    nsource = 1;
end
if nargin<4||isempty(type)
    type = 'ug';
end
nd = nsource; % number of components per diople
Covy = (Covy+Covy')/2;
% Covy = Covy*mean(csvd(L*L'))*nchan/trace(Covy);
[U,S,~] = svd(Covy);

S = diag(S);
rankn = rank(Covy);
% rankn = sum(S./S(1)>10*nchan*eps);
U = U(:,1:rankn);
S = S(1:rankn)+0.05*mean(S);
iCov = U*diag(1./S)*U';
% iCov = inv(Covy+0.1*trace(Covy)/nchan*eye(nchan));
% iCov = pinv(Covy);
Wimaging = zeros(size(L'));
switch type
    case 'ug'
        for iter = 1:nd:size(L,2)
            curIdx = iter:iter+nd-1;
            ft = pinv(L(:,curIdx)'*iCov*L(:,curIdx))*L(:,curIdx)'*iCov;
            Wimaging(curIdx,:) = ft; 
        end
        
    case 'ua'   
        for iter = 1:nd:size(L,2)
            curIdx = iter:iter+nd-1;
            fnt = L(:,curIdx)./norm(L(:,curIdx));
            ft = pinv(fnt'*iCov*fnt)*fnt'*iCov;
            Wimaging(curIdx,:) = ft; 
        end
        
    case 'ung'
        for iter = 1:nd:size(L,2)
            curIdx = iter:iter+nd-1;
            fnt = pinv(L(:,curIdx)'*iCov*L(:,curIdx));
            gamma = fnt*(L(:,curIdx)'*iCov^2*L(:,curIdx))*fnt;            
            Wimaging(curIdx,:) = diag(1./sqrt(diag(gamma)))*fnt*L(:,curIdx)'*iCov; 
        end  
end
end