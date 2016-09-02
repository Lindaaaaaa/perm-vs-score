% This function implements the score test for GxE interaction
% Input:    y is an nSub x 1 vector for phenotypes
%           X is an nSub x Ncov matrix for covariates
%           K is an nSub x nSub x nKer kernel matrix 
%           tol is the tolerance for the convergence of the ReML algorithm
%           maxIter is the maximum number of iterations for the ReML algorithm
%           alg = 1:average information, 2:expected information, 3:observed information
% Output:   flag = 1:ReML converged, 0:ReML did not converged
%           pval is the p-value
function [flag, pval] = SKATPermTestmultiKernel(y, X, K, nPerm,tol, maxIter, alg)


% Set parameters
nSub = length(y);
nKer = size(K,3);
K(:,:,end+1) = eye(nSub);

% Test each tau_i
pval = zeros(nKer,1); 
for i = 1:nKer 
    ind = find(1:nKer~=i);
    
    % Fit the null model using ReML algorithm
    [flag(i),tauSq] = REML(y, X, K(:,:,ind), tol, maxIter, alg); % Important to remove Ki when testing tau_i
    
    % Estimate covariance
    V = zeros(nSub,nSub);
    for j = 1:nKer
        if j < nKer % tau_i2
            V = V+tauSq(j)*K(:,:,ind(j));
        else % sigma2
            V = V+tauSq(j)*K(:,:,end);
        end
    end
    
    % Estimate projection matrix
    P = (eye(nSub)-(V\X)*((X'*(V\X))\X'))/V;

    % Estimate test statistic
    stat = y'*P*K(:,:,i)*P*y/2;
    
    % inference
    statPerm =nan(nPerm,1);
    for k=1:nPerm
        yperm=y(randperm(nSub));
        statPerm(k)=yperm'*P*K(:,:,i)*P*yperm/2;
    end
    pval(i) = (sum(statPerm>stat)+1)/(nPerm+1);
end

end
