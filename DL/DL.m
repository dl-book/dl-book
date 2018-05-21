% Copyright (c) 2016-2018 Paul Irofti <paul@irofti.net>
% 
% Permission to use, copy, modify, and/or distribute this software for any
% purpose with or without fee is hereby granted, provided that the above
% copyright notice and this permission notice appear in all copies.
% 
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

function [Dout, Xout, err, errextra, shared] = ...
    DL(Y, D, s, iternum, varargin)
%% Dictionary Learning (DL) Framework
% INPUTS:
%	Y -- signal set
%   D -- dictionary
%   s -- sparity target
%   iternum -- DL iterations
%
% OPTIONAL INPUTS:
%   upfunc -- dictionary update function pointer (default: K-SVD).
%             The function has to respect the API: 
%                [D,X] = func(Y,D,X,iter,replatoms,varargin)
%             Look in the dict/ directory for existing implementations.
%   upopts -- the update function specific options (passed on to upfunc)
%   spfunc -- sparse representation function pointer (default: OMP).
%             The function has to respect the API: 
%                X = func(Y,D,s,varargin)
%             Look in the reps/ directory for existing implementations.
%   spopts -- the representation function specific options 
%             (passed on to spfunc)
%
% PARAMETERS (key-value pairs):
%   replatoms -- unused dictionary atoms replacement strategy 
%               (default: 'zero')
%   postopts -- extra steps to perform after a DL iteration
%   postoptsargs -- arguments for extra DL iter steps
%   erropts -- compute other error reduction criteria than RMSE
%   outopts -- decide which dictionary and representations to return
%   shared -- shared data between sparse coding and dictionary update
%
% OUTPUTS:
%   Dout -- the resulting dictionary
%   Xout -- the corresponding sparse representations
%   err -- the RMSE at each iteration
%   errextra -- the custom error calculated at each iter (see erropts)
    
    p = inputParser;
    
    p.addOptional('upfunc', @ksvd);
    p.addOptional('upopts', {});
    p.addParameter('replatoms', 'zero', @(x) true);
    p.addParameter('spfunc', @omp);
    p.addParameter('spopts', {}, @(x) iscell(x));
    p.addParameter('initopts', {}, @(x) iscell(x));
    p.addParameter('postopts', '', @(x) true);
    p.addParameter('postoptsargs', {}, @(x) iscell(x));
    p.addParameter('erropts', {}, @(x) iscell(x));
    p.addParameter('outopts', 'best', @(x) true);
    p.addParameter('data_mask', []);  % !!! BD
    p.addParameter('shared', {}, @(x) iscell(x));

    p.parse(varargin{:});
    res = p.Results;
    shared = res.shared;

    Ymask = res.data_mask;  % !!! BD
    if ~isempty(Ymask)
      res.spopts = {Ymask};
      res.upopts = {Ymask};
      % res.erropts = {res.erropts{:}, 'data_mask', Ymask};
      res.erropts = {'', 'data_mask', Ymask};
    end

    bestrmse = Inf;
    err = zeros(1,iternum);
    errextra = zeros(1,iternum);
    
    [Y, D, X] = DL_init(Y, D, s, iternum, res.initopts);

    for iter = 1:iternum
        if strcmp(func2str(res.spfunc), 'NOP') == 0
            [X, shared] = res.spfunc(Y, D, s, shared, res.spopts{:});
        end
        [D, X, shared] = res.upfunc(Y, D, X, iter, ...
            res.replatoms, shared, res.upopts{:});
                
        if ~isempty(res.postopts)
            [D,X] = DL_extra(Y, D, X, iter, iternum, ...
                res.postopts, res.postoptsargs{:});
        end
        [err(iter), errextra(iter)] = ...
            DL_error(Y, D, X, iter, shared, res.erropts{:});
        
        % Store the best dictionary
        switch res.outopts
            case 'best'
                if (err(iter) < bestrmse)
                    bestrmse = err(iter);
                    Dout = D;
                    Xout = X;
                end
            case 'last'
                Dout = D;
                Xout = X;
        end
    end
end

function [Y, D, X] = DL_init(Y, D, s, iternum, initopts)
    p = inputParser();
    
    p.addParameter('initreps', 0);
    
	p.parse(initopts{:});
    
    X = p.Results.initreps;
    
    if iscell(D)
        N = size(Y,2);
        p1 = size(D{1},1);
        p2 = size(D{2},2);
        Y = reshape(Y,p1,p2,N);
    end
end

function [D,X] = DL_extra(Y, D, X, iter, iternum, postops, varargin)
    persistent F;
    persistent mutual_coh_ipr;
    persistent nit_ipr;

    if iter == 1
        p = inputParser();
        p.KeepUnmatched=true;
        p.addParameter('mutual_coh_ipr', 0.2);
        p.addParameter('nit_ipr', 5);
        p.parse(varargin{:});
        mutual_coh_ipr = p.Results.mutual_coh_ipr;
        nit_ipr = p.Results.nit_ipr;
    end

    switch postops
        case 'worstrep'
            r = norms(Y-D*X);
            ndx = sum(abs(X),2) == 0;
            na = sum(ndx);
            if na > 0
                [~,ndx_sort] = sort(r, 'descend');
                D_new = Y(:,ndx_sort(1:na));
                D_new = D_new./repmat(norms(D_new), size(D,1), 1);
                D(:,ndx) = D_new;
            end
        case 'atomsusage'
            if isempty(F)
                F = getframe;
            end
            F(iter) = atomsusage(D,X);
            if iter == iternum
                movie(F, 20);
            end
        case 'finalLS'
            if iter == iternum
                X = ompreg(Y,D,X,0);
            end
        case 'ipr'   % make incoherent dictionary
            D = ipr(Y, D, X, mutual_coh_ipr, nit_ipr);
    end
end

function [rmse, errextra] = DL_error(Y, D, X, iter, shared, varargin)
    persistent coh;
    persistent vanish;
    persistent reg;
    persistent erropts;
    persistent Q2s;
    persistent Dtrue;
    %persistent Ymean;
    persistent Ymask;

    if iter == 1
        p = inputParser();
        p.KeepUnmatched=true;
        p.addOptional('erropts', '', @(x) true);
        p.addParameter('coh', 6);
        p.addParameter('reg', 0.01);
        p.addParameter('vanish', 0.95);
        p.addParameter('Q2s', eye(size(Y,1)));
        p.addParameter('Dtrue', zeros(size(D)));
        %p.addParameter('Ymean', zeros(size(Y)));
        p.addParameter('data_mask', []);  % !!! BD
        p.parse(varargin{:});
        erropts = p.Results.erropts;
        coh = p.Results.coh;
        reg = p.Results.reg;
        vanish = p.Results.vanish;
        Q2s = p.Results.Q2s;
        Dtrue = p.Results.Dtrue;
        %Ymean = repmat(p.Results.Ymean,size(Y,1),1);
        Ymask = p.Results.data_mask;
    end

    if ndims(Y) == 3
        [p1, p2, N] = size(Y);
        Nerr = zeros(N,1);
        for k = 1:N
            Nerr(k) = norm(Y(:,:,k) - D{1}*X(:,:,k)*D{2}, 'fro')^2;
        end
        rmse = sqrt(sum(Nerr))/sqrt(p1*p2*N);
        errextra = 0;
        return;
    end

    [p, m] = size(Y);
    
    if ~isempty(Ymask)
      rmse = norm(Ymask .* (Y - D*X), 'fro') / sqrt(sum(sum(Ymask)));  % !!! BD
      errextra = norm(Y - D*X, 'fro') / sqrt(p*m); % !!! BD
      return;
    end

    switch erropts
        case 'coherr' 
            n = size(D,2);
            errextra = norm(D*X-Y,'fro')^2+coh*norm(D'*D-eye(n),'fro')^2;
        case 'regerr'
            reg = reg * vanish;
            errextra = norm(D*X-Y,'fro')^2 + reg*norm(X,'fro')^2;
        case '2s'
            rmse = norm(Y - Q2s*D*X, 'fro') / sqrt(p*m);
            errextra = 0;
            return;
        case 'var2s'
            Q = shared{1};
            rmse = norm(Y - Q*D*X, 'fro') / sqrt(p*m);
            errextra = 0;
            return;
        case 'avgnnz'
            errextra = mean(sum(X ~= 0, 1));
        case 'recov'
            errextra = recovered_atoms(D, Dtrue);
        otherwise
            errextra = 0;
    end
    
    %% Errror
    rmse = norm(Y - D*X, 'fro') / sqrt(p*m);
    %if strcmp(erropts, 'relrmse')
    %    rmse = norm(Y+Ymean - (D*X+Ymean), 'fro') / sqrt(p*m);
    %end

end
