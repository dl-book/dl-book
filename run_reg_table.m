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

%% Table 4.1 data: RMSE for standard and regularized DL
clear; clc; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
m = 64;              % problem dimension
n = 256;             % number of atoms in the dictionary
N = 4096;            % number of training signals
ss = [4 6 8 10 12];  % sparsity constraint
reg = 0.05;          % regularization
vanish = 1;          % regularization vanishing factor
pa = n;              % number of atoms to update in parallel
iters = 50;          % DL iterations
rounds = 10;         % test rounds
itN = 1;             % SimCO internal iterations

generate_data = true;   % generate initial data
ots = '20170818164513'; % else use data timestamp

% Dictionary update routines
updates = {'ksvd', 'ksvd_reg', 'aksvd', 'aksvd_reg', 'simco', 'simco_reg'};
% Unused atoms replacement strategy
replatom = 'no';
%%-------------------------------------------------------------------------
datadir = 'data\';
dataprefix = 'tab_4_1_reg';
imdir = 'img\';
images = {'barbara.png', 'boat.png', 'house.png', 'lena.png', 'peppers.png'};
%%-------------------------------------------------------------------------
timestamp = datestr(now, 'yyyymmddHHMMss');
if generate_data
    Y = cell(rounds,1);
    D0 = cell(rounds,1);
    for r = 1:rounds
        Y{r} = mkimgsigs(imdir,images,m,N,'distinct');
        D0{r} = normc(randn(m,n));
    end
    save([datadir 'tab_4_1_reg_init-' timestamp '.mat'], 'Y', 'D0');
else
    load([datadir 'tab_4_1_reg_init-' otime '.mat'], 'Y', 'D0');
end


methods = length(updates);
    
params = {'reg', reg, 'vanish', vanish, 'itN', itN};
erropts = {'regerr', 'reg', reg};

for s = ss
    fprintf('s=%d\n', s); 
    
    %% Rounds
    Yr = zeros(m,N);
    D0r = zeros(m,n);
    Dall = zeros(rounds,methods, m, n);
    Xall = zeros(rounds,methods, n, N);
    errs = zeros(rounds,methods, iters);
    criteria = zeros(rounds,methods, iters);
    
    for r = 1:rounds
        fprintf('%d', mod(r, 10)); 
        
        D0r(:,:) = D0{r};
        Yr(:,:) = Y{r};

        for j = 1:methods
            [Dall(r,j,:,:), Xall(r,j,:,:), errs(r,j,:), criteria(r,j,:)] = ...
                DL(Yr, D0r, s, iters, str2func(updates{j}), params, ...
                'replatom', replatom, 'erropts', erropts);
        end
    end
    %% Write out data
    matfile = sprintf('%s%s-m%d-n%d-N%d-s%d-i%d-%s.mat', ...
         datadir, dataprefix, m, n, N, s, iters, timestamp);    
    save(matfile, 'Dall', 'Xall', 'errs', 'criteria');
end % sparsity loop