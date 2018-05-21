function [X, shared] = omp_ker(K, Kz, A, s, shared, varargin)
% Kernel Orthogonal Matching Pursuit algorithm
% Can be used inside a DL algorithm or on its own
% Input:
%   K       - kernel matrix
%   Kz      - kernel matrix between training and test signals (in DL, Kz==K)
%             if empty, then Kz = K
%   A       - current dictionary
%   s       - sparsity constraint
%
% Output:
%   X       - sparse representations

% BD 22.12.2017

    if isempty(varargin)
        ompparams = {'checkdict', 'off'};
        do_sparse = true;
    else
        if strcmp(varargin{1}, 'error')
            do_sparse = false;
            error = varargin{2};
            if length(varargin) > 2
                ompparams = varargin{3:end};
            else
                ompparams = {'checkdict', 'off'};
            end
        else
            do_sparse = true;
            ompparams = varargin;
        end
    end
    if isempty(Kz)      % when OMP applied on training signals
      Kz = K;
    end
    if do_sparse
        X = omp_sparse(A'*Kz, A'*K*A, s, ompparams{:});
    else
        X = omp_err(A'*Kz, sum(K), A'*K*A, error, ompparams{:});
    end
end
