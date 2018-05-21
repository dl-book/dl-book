% Copyright (c) 2017-2018 Paul Irofti <paul@irofti.net>
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

%% Table 4.1 extra SimCO itN data
clear; clc; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
m = 64;                       % problem dimension
n = 256;                      % number of atoms in the dictionary
N = 4096;                     % number of training signals
ss = [4 6 8 10 12];           % sparsity constraint
reg = 0.05;                   % regularization
vanish = 1;                   % regularization vanishing factor
iters = 50;                   % DL iterations
rounds = 10;                  % test rounds
iitN = [10, 20, 30, 40, 50];  % SimCO internal iterations

ots = '20170818164513'; % timestamp, copy the one from run_reg_table

% Dictionary update routines
updates = {'simco_reg'};
% Unused atoms replacement strategy
replatom = 'no'; 
%%-------------------------------------------------------------------------
datadir = 'data\';
dataprefix = 'tab_4_1_reg';
%%-------------------------------------------------------------------------
load([datadir 'tab_4_1_reg_init-' ots '.mat'], 'Y', 'D0');

methods = length(updates);

for itN = iitN
    params = {'reg', reg, 'vanish', vanish, 'itN', itN};
    erropts = {'regerr', 'reg', reg};
    fprintf('itN=%d\n', itN);

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
        matfile = sprintf('%s%s-itN%d-m%d-n%d-N%d-s%d-i%d-%s.mat', ...
             datadir,dataprefix, itN, m, n, N, s, iters, ots);
        save(matfile, 'Dall', 'Xall', 'errs', 'criteria');
    end % sparsity loop
end % itN