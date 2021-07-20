function Wspatial=CommonSpatialPattern(X1,X2,type)
% common spatial pattern
% Input 
%   X1, X2:EEG data, Channel*Sample*Trial
% Output
%   spatial filter
if nargin<3
    type = 1;
end
nchan = size(X1,1);
nsam = size(X1,2);
numX1 = size(X1,3);
numX2 = size(X2,3);

covX1 = zeros(size(X1,1),size(X1,1),size(X1,3));
covX2 = zeros(size(X2,1),size(X2,1),size(X2,3));
for num = 1:numX1
    temp = squeeze(X1(:,:,num));
    temp = bsxfun(@minus,temp,mean(temp,2));
    temp = temp*temp';
    covX1(:,:,num) = temp/trace(temp);
end
for num = 1:numX2
    temp = squeeze(X2(:,:,num));
    temp = bsxfun(@minus,temp,mean(temp,2));
    temp = temp*temp';
    covX2(:,:,num) = temp/trace(temp);
end
covX1 = mean(covX1,3);
covX2 = mean(covX2,3);

switch type
    case 2
        c1 = 1-2/nchan;
        c2 = nsam+1-2/nchan;
        c3 = 1-nsam/nchan;        
        rho1 = (c1*trace(covX1^2) + trace(covX1)^2) / (c2*trace(covX1^2) + c3*trace(covX1)^2);
        covX1=(1-rho1)*covX1+rho1*trace(covX1)*eye(nchan);
        rho2 = (c1*trace(covX2^2) + trace(covX2)^2) / (c2*trace(covX2^2) + c3*trace(covX2)^2);
        covX2=(1-rho2)*covX2+rho2*trace(covX2)*eye(nchan);
        [U,S,~] = eig(covX1,covX2);
        [~,idx] = sort(diag(S),'descend');
        Wspatial = U(:,idx);    
    case 3
        Wspatial = zeros(nchan);
        alpha = 1e-3*trace(covX2);
        [U,S,~] = eig(covX1,covX2+alpha*eye(nchan));
        [~,idx] = sort(diag(S),'descend');
        U = U(:,idx);
        Wspatial(:,1:floor(nchan/2)) = U(:,1:floor(nchan/2));
        [U,S,~] = eig(covX2,covX1+alpha*eye(nchan));
        [~,idx] = sort(diag(S));
        U = U(:,idx);  
        Wspatial(:,floor(nchan/2)+1:end) = U(:,floor(nchan/2)+1:end);
    otherwise
        [U,S,~] = eig(covX1,covX2);
        [~,idx] = sort(diag(S),'descend');
        Wspatial = U(:,idx);
end

for i = 1:size(Wspatial,2)
    Wspatial(:,i) = Wspatial(:,i)./norm(Wspatial(:,i));
end

end



%%% shrinkage
% OAS describer in:
%     [1] "Shrinkage Algorithms for MMSE Covariance Estimation"
%     Chen et al., IEEE Trans. on Sign. Proc., Volume 58, Issue 10, October 2010.
% nchan = size(X1,1);
% nsam = size(X1,2);
% c1 = 1-2/nchan;
% c2 = nsam+1-2/nchan;
% c3 = 1-nsam/nchan;
% 
% rho1 = (c1*trace(covX1^2) + trace(covX1)^2) / (c2*trace(covX1^2) + c3*trace(covX1)^2);
% covX1=(1-rho1)*covX1+rho1*trace(covX1)*eye(nchan);
% rho2 = (c1*trace(covX2^2) + trace(covX2)^2) / (c2*trace(covX2^2) + c3*trace(covX2)^2);
% covX2=(1-rho2)*covX2+rho2*trace(covX2)*eye(nchan);

% TRCSP
%   "Regularizing Common Spatial Patterns to Improve BCI Designs:Unified Theory and New Algorithms"
