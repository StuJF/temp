function Wimaging = MNFilt(L,lamda,Covy,nsource,type,varargin)
% W = MNFilt(L,lamda,Covy,type,nsource)
% Minimum-norm method for EEG Source Imaging
% Input:
%   L: Gain matrix
%   Covy: data covariance (default is lamda*L*L'+I)
%   lamda: regularization parameters (default is 1/SNR^2)
%   type: min-norm method (default is sloreta)

[nchan,~] = size(L);
remat = eye(nchan) - ones(nchan)/nchan;
if nargin<2||isempty(lamda)
    lamda = 0.1;
end

if nargin<3||isempty(Covy)    
    Covy = remat;
end
if nargin<4||isempty(nsource)
    nsource = 1;
end
if nargin<5||isempty(type)
    type = 'sloreta';
end

Covs = 1; % prior covariance of source
nd = nsource; % number of components per diople
%Covy = diag(Covy);
Covy = (Covy+Covy')/2;
[U,S,~] = svd(Covy);
Sn = sqrt(diag(S));
Sn = Sn(1:rank(Covy));
U = U(:,1:rank(Covy));
iW = U*diag(1./Sn)*U'; % whitener
if ~(strcmp(type,'sloreta')||strcmp(type,'eloreta'))
    L = WeightGain(L,nd);
end
%iW = eye(nchan);
L = iW*L;
[UL,SL,~] = svd(L*Covs*L');
SL = diag(SL);
%lamda = lamda*mean(SL);
iG = UL*diag(1./(SL+lamda*mean(SL)))*UL';
%iG = inv(L*Covs*L'+lamda*eye(nchan));
switch type
    case 'mne'
        Wimaging = Covs*L'*iG*iW;
    case 'sloreta'
        Wimaging = sloreta_filt(L,Covs,iG,nd)*iW;
    case 'eloreta'
        Wimaging = eloreta_filt(L,Covs,iG,lamda,nd,30)*iW;
    case 'dspm'
        Kernel = Covs*L'*iG*iW;
        R = sum(Kernel.^2,2);
        if nd>1
            R = sum(reshape(R,nd,[]),1);
            R =  repmat(R,[nd,1]);
            R = R(:);
        end
        Wimaging = bsxfun(@rdivide,Kernel,sqrt(R));
        
end

end

function L = WeightGain(L,nd)
normL = sqrt(sum(L.^2));
if nd==1
    Lbound = max(normL)/100;
    W = normL;
    W(W<Lbound) = Lbound;
    W = max(normL)./W;
    L = bsxfun(@times,L,sqrt(W));
end
if nd>1
    normL = reshape(normL,3,[]);
    normL = sum(normL,1);
    Lbound = max(normL)/100;
    W = normL;
    W(W<Lbound) = Lbound;
    W = max(normL)./W;
    W = repmat(W,nd,1);
    W = reshape(W,1,[]);
    L = bsxfun(@times,L,sqrt(W));
end

end

function W = sloreta_filt(L,Covs,iG,nd)
W = Covs*L'*iG;
if nd==1   
    W = bsxfun(@rdivide,W,sqrt(sum(W.*L',2)));
    return
end
if nd>1
    for iter = 1:nd:size(W,1)
        curIdx = iter:iter+nd-1;
        R = W(curIdx,:)*L(:,curIdx);
        W(curIdx,:) = sqrtm(pinv(R))*W(curIdx,:);
        %W(curIdx,:) = diag(1./sqrt(diag(R)))*W(curIdx,:);
    end
    return
end
end

function Wout = eloreta_filt(L,Covs,iG,lamda,nd,maxIter)
nchan = size(L,1);
nc = size(L,2)/nd;
Wold = repmat(eye(nd),1,nc);
Wold = reshape(Wold,nd,nd,[]);
W = Wold;
Winv = W;
M = iG;
for iter = 1:maxIter
    for ic = 1:nc
        curIdx = (1+(ic-1)*nd):ic*nd;
        W(:,:,ic) = sqrtm(L(:,curIdx)'*M*L(:,curIdx));  
        Winv(:,:,ic) = pinv(W(:,:,ic));
    end
    winvk = zeros(size(L'));
    for ic=1:nc
       curIdx = (1+(ic-1)*nd):ic*nd;
       winvk(curIdx,:) = Winv(:,:,ic)*L(:,curIdx)';
    end
    kwinvk = L*winvk;
    M = pinv(kwinvk+lamda*trace(kwinvk)*eye(nchan)/nchan); 
    if norm(reshape(W-Wold,[],1))/norm(reshape(Wold,[],1))<1e-6
        break;
    end
    Wold = W;   
end
Wout = winvk*M;
% nchan = size(L,1);
% %Wold = eye(size(L,2));
% %Wold = reshape;
% M = iG;
% if nd==1
%     Wold = ones(size(L,2),1);
%     M = sqrtm(M);
%     for iter = 1:maxIter
%         W = sqrt(sum((L'*M).^2,2));
%         M = sqrtm(pinv(bsxfun(@rdivide,L,W)*L'+lamda*eye(nchan)));
%         if norm(W-Wold)/norm(Wold)<1e-6
%             break;
%         end
%         Wold = W;
%     end
%     Wout = bsxfun(@rdivide,L'*M,W);
%     return;
% end
% 
% if nd>1
%     Wold = repmat(eye(nd),size(L,2)/nd,1);
%     W = Wold;
%     for iter = 1:maxIter
%         Winv = zeros(size(L,1));
%         for num = 1:nd:size(L,2)
%             curIdx = num:num+nd-1;
%             W(curIdx,:) = L(:,curIdx)'*M*L(:,curIdx);
%             W(curIdx,:) = sqrtm(W(curIdx,:));
%             
%             % trace(W(curIdx,curIdx))*eye(nd)/10^10);
%         end
%         
%         for num = 1:nd:size(L,2)
%             curIdx = num:num+nd-1;
%             Winv = Winv+L(:,curIdx)*pinv(W(curIdx,:))*L(:,curIdx)';
%         end
%         M = pinv(Winv+lamda*eye(nchan));
%         if norm(W-Wold)/norm(Wold)<1e-6
%             break;
%         end
%         Wold = W;
%     end
%     Wout = zeros(size(L'));
%     for num = 1:nd:size(L,2)
%         curIdx = num:num+nd-1;
%         Wout(curIdx,:) = pinv(W(curIdx,:))*L(:,curIdx)';
%     end
%     Wout = Wout*M;
% end


end