function [SIdx,model,AIC]=LARS(X,Y,K)
% least angle regression
% Input
%   X:features with 0 mean and unit length;each row belongs to one
%       class
%   Y:label with 0 mean
%   K:number of feature to select

%  Output
%   SIdx:index of selected features
%   model:weights for X i.e. for ith iteration result Ypred = X*model(i,:)'
%   AIC:Akaike information criterion
if nargin<3
    K=size(X,2);
end
delta = 1e-4;

X = X-mean(X,1);
nm = sqrt(sum(X.^2));
X = bsxfun(@rdivide,X,nm); % normalize
X(:,nm==0)=0;
Y = Y-mean(Y); % centering

yC = zeros(size(Y)); % current estimation of y
SIdx = []; 
AIC = zeros(1,K); %AIC (Akaike information criterion)
weights = zeros(K,1);
model = [];
for iter = 1:K
   Coerr = X'*(Y - yC); % current correlation
   rIdx = setdiff(1:K,SIdx); % index that not belong to SIdx
   [Cmax,idxC] = max(abs(Coerr(rIdx))); % find max correlation
   if norm(Y-yC)<=delta
       break;
   end
   SIdx = [SIdx rIdx(idxC)]; % add new index
   signC = sign(Coerr(SIdx)); % sign of correlation
   XA = X(:,SIdx)*diag(signC);
   G = pinv(XA'*XA);
   oneA = ones(length(SIdx),1);
   Aa = (oneA'*G*oneA)^(-0.5);
   wA = Aa*G*oneA;
   uA = XA*wA;
   a = X'*uA;
   if iter < length(Y)
       M = [(Cmax-Coerr(SIdx))./(Aa-a(SIdx));(Cmax+Coerr(SIdx))./(Aa+a(SIdx))];
       M(M<=0) = +inf;
       gamma = min(M);
   else
       gamma = Cmax/Aa;
   end
   yC = yC+gamma*uA;
   weights(SIdx) = weights(SIdx)+gamma*diag(signC)*wA; 
   model = [model;weights'];
   AIC(iter) = 2*iter+size(X,1)* log(sumsqr(Y-X*weights)/size(X,1));
   
end
    
