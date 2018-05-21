function [f,X,g,freal] = fg_tilde_eval01(Y,D,Omega,IPara)

% fg_tilde_eval01 computes the gradient descent direction for LineSearch 
% in regularized SimCO version
%
% References:
% W. Dai, T. Xu, and W. Wang, 
% "Simultaneous Codeword Optimization (SimCO) for Dictionary Update and Learning,"
% submitted to IEEE Transactions on Signal Processing, October 2011.
% Full text is available at http://arxiv.org/abs/1109.5302
%
%
% Wei Dai, Tao Xu, Wenwu Wang
% Imperial College London, University of Surrey
% wei.dai1@imperial.ac.uk  t.xu@surrey.ac.uk  w.wang@surrey.ac.uk
% October 2011

[m,n] = size(Y);
d = size(D,2);
X = zeros(d,n);
OmegaL = sum(Omega,1);
mu = IPara.mu; %the parameter of regularized item
mu_sqrt = sqrt(mu);
for cn = 1:n
    L = OmegaL(cn);
    X(Omega(:,cn),cn) = ...
        [D(:,Omega(:,cn));diag(mu_sqrt*ones(L,1))] \ [Y(:,cn);zeros(L,1)];
end
Yr = Y - D*X;
% the cost function with regularized term
f = sum(sum(Yr.*Yr))+mu*sum(sum(X.*X));
freal = sum(sum(Yr.*Yr));
if nargout >= 3
    g = -2*Yr*X';
    % additional steps to make sure the orthoganilty
    DGcorr = sum(D.*g,1);
    g = g - D.*repmat(DGcorr,m,1);
end
