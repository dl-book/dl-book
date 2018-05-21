function [D,X,OPara] = DictUpdate03 (Y, D, X, IPara)
% DictUpdate03 is the dictionary update function in SimCO.
% Given the initial dictionary D, initial sparse coefficient matrix X and
% the traning data matrix Y, this function produces the updated D and X
% through itN iterations of line search algorithm in DictLineSearch03
% 
% References:
% W. Dai, T. Xu, and W. Wang, 
% "Simultaneous Codeword Optimization (SimCO) for Dictionary Update and Learning,"
% submitted to IEEE Transactions on Signal Processing, October 2011.
% Full text is available at http://arxiv.org/abs/1109.5302
%
% See also K_DictUpdate02
%
% Wei Dai, Tao Xu, Wenwu Wang
% Imperial College London, University of Surrey
% wei.dai1@imperial.ac.uk  t.xu@surrey.ac.uk  w.wang@surrey.ac.uk
% October 2011



% IPara.I : the subset of training atoms in SimCO. Default value is the universal set.
% IPara.dispN : the display interval. Default value is 20.
% IPara.DebugFlag : the flag for debugging ill-conditioned problem.
%   default value is 0. Set it to 1 for debugging 
% IPara.itN : the number of iterations. default value is 100 for each
%   values of mu(the parameter of regularized item) in synthetic tests
% IPara.gmin : if the gradient frobenius norm is less than gmin, then the
%   gradient can be viewed as zero. default value is 1e-5.
% IPara.Lmin : t4-t1 should be larger than Lmin. Default value is 1e-6.
% IPara.rNmax : the number of iterative refinement in Part B. default value 
%   is 3.
% IPara.t4 : the initial value of t4. Default value is 1e-2.
% IPara.mu : regularization parameter, if mu != 0, regularized SimCO, otherwise prime SimCO. 

% If D = [], we generate D randomly

% OPara.f0 : initial value of f
% OPara.f1 : final value of f
% OPara.f0real : true initial value of cost function without the effect of mu
% OPara.f1real : true final value of cost function without the effect of mu
% OPara.topt : optimal t value
% OPara.gn2 : the frobenius norm of the gradient
% OPara.Flag = 1 : gradient is too small
% OPara.Flag = 0 : successful update
% OPara.Flag = -1 : t is too small
% OPara.CondNum : the condition number of trained dictionary for each
%   iteration


%% initialization

I = IPara.I;
DebugFlag = IPara.DebugFlag; 
itN = IPara.itN;  
OPara.Flag = zeros(1,itN);
OPara.f0 = zeros(1,itN);
OPara.f1 = zeros(1,itN);
OPara.f0real = zeros(1,itN);
OPara.f1real = zeros(1,itN);
OPara.gn2 = zeros(1,itN);
OPara.topt = zeros(1,itN);
d = size(X,1);
m = size(Y,1);

% debug the code with condition number of the trained dictionary
if DebugFlag == 1
    OPara.CondNum = zeros(1,itN);
    OPara.CondNum = OPara.CondNum(1:itN);
end

if isempty(D)
    D = randn(m,d);
    Dn2 = sqrt(sum(D.*D,1));
    D = D.*repmat(1./Dn2,m,1);
end

%Ic is the complementary set of I
I = intersect(1:d,I);
Ic = setdiff(1:d,I);
Yp = Y - D(:,Ic)*X(Ic,:);
Omega = X~=0;
ColUpdate = sum(Omega(I,:),1)~=0;
YI = Yp(:,ColUpdate);
DI = D(:,I);
XI = X(I,ColUpdate);
OmegaI = Omega(I,ColUpdate);
f_YIComp = norm( Yp(:,~ColUpdate),'fro' )^2;

%% gradient descent line search
for itn = 1:itN
    if itn == 1
        OPara.f0real(itn) = norm(Y-D*X,'fro')^2;
    else
        OPara.f0real(itn) = OPara.f1real(itn-1);
    end
    
    %us the line search mechanism for dictionary update
    [DI,XI,OParaLS] = DictLineSearch03(YI,DI,OmegaI,IPara);
    D(:,I) = DI;
    X(I,ColUpdate) = XI;    
    OPara.Flag(itn) = OParaLS.Flag;
    OPara.f0(itn) = OParaLS.f0+f_YIComp;
    OPara.f1(itn) = OParaLS.f1+f_YIComp;
    OPara.f1real(itn) = norm(Y-D*X,'fro')^2;
    OPara.gn2(itn) = OParaLS.gn2;
    OPara.topt(itn) = OParaLS.topt;
    if DebugFlag == 1
        D(:,I) = DI;
        Dc = Dcondition(Y,D,Omega);
        OPara.CondNum(itn) = max(Dc(1,:));
    end
    
    if OParaLS.Flag ~= 0
        OPara.Flag = OPara.Flag(1:(itn));
        OPara.f0 = OPara.f0(1:(itn));
        OPara.f1 = OPara.f1(1:(itn));            
        OPara.f0real = OPara.f0real(1:(itn));
        OPara.f1real = OPara.f1real(1:(itn));
        OPara.gn2 = OPara.gn2(1:(itn));
        OPara.topt = OPara.topt(1:(itn));
        if DebugFlag == 1
            OPara.CondNum = OPara.CondNum(1:(itn));
        end  
        break;
    end
    IPara.t4 = OParaLS.topt;
    
end

% finalize
D(:,I) = DI;
X(I,ColUpdate) = XI;
