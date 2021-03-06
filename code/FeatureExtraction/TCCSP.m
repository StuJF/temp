function [Wspatial,Lap] = TCCSP(X1,X2,type,tao)
% modified temperal constrained common spatial pattern
% Input 
%   X1, X2:EEG data, Channel*Sample*Trial
% Output
%   spatial filter, laplacian matrix
if nargin<3
    type = 1;
end
if nargin<4
    tao = 1;
end
nchan = size(X1,1);
nsam = size(X1,2);
numX1 = size(X1,3);
numX2 = size(X2,3);

mean_X1 = zeros(nchan,nsam);
mean_X2 = zeros(nchan,nsam);
for num = 1:numX1
    temp = squeeze(X1(:,:,num));
    temp = bsxfun(@minus,temp,mean(temp,2));
    temp = temp/sqrt(trace(temp*temp'));
    mean_X1 = mean_X1+temp;
end
mean_X1 = mean_X1/num;
for num = 1:numX2
    temp = squeeze(X2(:,:,num));
    temp = bsxfun(@minus,temp,mean(temp,2));
    temp = temp/sqrt(trace(temp*temp'));
    mean_X2 = mean_X2+temp;
end
mean_X2 = mean_X2/num;

W_X1 = sparse(zeros(nsam));
W_X2 = W_X1;
switch type
    case 2
        for itp = 1:nsam
           for jtp = max(1,itp-tao):min(itp+tao,nsam) 
                W_X1(itp,jtp) = exp(-norm(mean_X1(:,itp)-mean_X1(:,jtp)));
                W_X2(itp,jtp) = exp(-norm(mean_X2(:,itp)-mean_X2(:,jtp)));
           end
        end
        
    otherwise
        for itp = 1:nsam
           for jtp = max(1,itp-tao):min(itp+tao,nsam) 
                W_X1(itp,jtp) = exp(corr(mean_X1(:,itp),mean_X1(:,jtp)));
                W_X2(itp,jtp) = exp(corr(mean_X2(:,itp),mean_X2(:,jtp)));
           end
        end        
end
L



end