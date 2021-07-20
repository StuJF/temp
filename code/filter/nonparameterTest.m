% 不太懂整个无参估计是用的什么方法
% fast calculation of kendall corr

for sub = 1:76
    disp(['The' num2str(sub) ' subject Positive symptom!']);
    Training_data = SCZ_Base;
    Training_scores = PANSS_base(:,1);
    
    % Select training data and testing data
    Testing_data = Training_data(:,:,sub);
    Testing_scores = Training_scores(sub);
    
    Training_data(:,:,sub) = [];
    Training_scores(sub) = [];
    
    Testing_data(isnan(Testing_data)) = 0;
    Training_data(isnan(Training_data)) = 0;
    
   
    permdata = permute(Training_data,[3 1 2]);
    nij = size(permdata,1);
    [rankMat,adjMat] = tiedrank(permdata,1);
    n2const = nij*(nij-1)/2;
    
    
    
    for N = 1:500
        RandIndex = randperm(76-1);
        Training_scores = Training_scores(RandIndex);
        [yrank,yadj] = tiedrank(Training_scores,1);
        K = zeros([size(permdata,2),size(permdata,3)]);
        for iter = 1:nij-1
            temp = rankMat(iter,:,:)-rankMat(iter+1:end,:,:);
            stemp = yrank(iter)-yrank(iter+1:end);
            K = K+squeeze(sum(sign(bsxfun(@times,temp,stemp)),1));
        end
        R = K./sqrt((n2const-yadj(1)).*squeeze(n2const-adjMat(1,:,:)));
        stdK = sqrt(n2const*(2*nij+5)./9)+adjMat(1,:,:).*yadj(1)./n2const ...
            + adjMat(2,:,:).*yadj(2)./(18*n2const*(nij-2))-(adjMat(3,:,:) + yadj(3))./18;
        muK = -(abs(K)-1)./stdK;
        P = erfc(-muK./sqrt(2));
        P(P>1) = 1;
%         for i = 1:246
%             for j = i+1:246
%                 [r, p] = corr(squeeze(Training_data(i,j,:)),Training_scores,'type','Kendall');
%                 R(i,j) = r; P(i,j) = p;
%             end
%         end
  
        P(isnan(P)) = 1; R(isnan(R)) = 0; 
        P_mask = (P.*(P<0.01)>0);
       
        r_po = (P_mask.*R)>0;
        r_ne = (P_mask.*R)<0;
    
        for s = 1:75
            S_po(s,1) = sum(sum(Training_data(:,:,s).*r_po))/2;        
            S_ne(s,1) = sum(sum(Training_data(:,:,s).*r_ne))/2;        
        end
    
        fit_pos = polyfit(S_po, Training_scores,1);
        fit_neg = polyfit(S_ne, Training_scores,1); 
    
        test_po = sum(sum(r_po.*Testing_data))/2;
        test_ne = sum(sum(r_ne.*Testing_data))/2;
        Predict_po_PANSS1(sub,N) = fit_pos(1)*test_po+fit_pos(2);
        Predict_ne_PANSS1(sub,N) = fit_neg(1)*test_ne+fit_neg(2);
    
        b = regress(Training_scores,[S_po,S_ne,ones(75,1)]);
        Predict_PANSS1(sub,N) = b(1)*test_po+b(2)*test_ne+b(3);
    
    end
end