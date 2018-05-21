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

%% Figure 4.1 data: RMSE for standard and regularized DL
clear; clc; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
m = 16;         % problem dimension
n = 32;         % number of atoms in the dictionary
NN = 64:32:512; % number of training signals
s = 4;          % sparsity constraint
reg = 0.1;      % regularization
vanish = 1;     % regularization vanishing factor
iters = 50;     % DL iterations
rounds = 50;    % test rounds
itN = 1;        % SimCO internal iterations

generate_data = true;   % generate initial data
ots = '20170818164513'; % else use data timestamp

datadir = 'data\';
dataprefix = 'fig_4_1_reg';

% Dictionary update routines
updates = {'ksvd', 'ksvd_reg', 'aksvd', 'aksvd_reg', 'simco', 'simco_reg'};
% Unused atoms replacement strategy
replatom = 'no';
%%-------------------------------------------------------------------------
timestamp = datestr(now, 'yyyymmddHHMMss');
methods = length(updates);
  
params = {'reg', reg, 'vanish', vanish, 'regstop', 31, 'itN', itN};
erropts = {'regerr', 'reg', reg};

if generate_data
    for N = NN
        D = zeros(rounds, m, n);
        Yawgn = zeros(rounds, m, N);
        for r = 1:rounds
            [D(r,:,:), ~, Yawgn(r,:,:)] = gen_synth_data(m,n,N,s,Inf);
        end
        matfile = sprintf('%s%s_init-m%d-n%d-N%d-s%d-i%d-%s.mat', ...
             datadir, dataprefix, m, n, N, s, iters, timestamp);
        save(matfile, 'D', 'Yawgn');
    end
    ots = timestamp;
end

for N = NN
    fprintf('m=%d\n', N); 
    
    matfile = sprintf('%s%s_init-m%d-n%d-N%d-s%d-i%d-%s.mat', ...
         datadir, dataprefix, m, n, N, s, iters, ots);
    load(matfile, 'Yawgn', 'D');
    
    %% Rounds
    Yr = zeros(m,N);
    D0r = zeros(m,n);
    Dall = zeros(rounds,methods, m, n);
    Xall = zeros(rounds,methods, n, N);
    errs = zeros(rounds,methods, iters);
    criteria = zeros(rounds,methods, iters);
    
    for r = 1:rounds
        fprintf('%d', mod(r, 10)); 
        
        D0r(:,:) = D(r, :, :);
        Yr(:,:) = Yawgn(r, :, :);

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
end