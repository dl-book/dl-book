% Copyright (c) 2017, 2018 Paul Irofti <paul@irofti.net>
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

%% Basic Algorithms: DL on images -- LAOLS tests
%% run_DL.m needs to be executed before this
clear; clc; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
m = 64;                 % problem dimension
nn = 128:64:512;        % number of atoms in the dictionary
N = 4000;               % number of training signals
ss = 6;                 % sparsity constraint
K = 100;                % DL iterations
rounds = 10;            % test rounds
Lahead = {3};           % look ahead
ots = '20170829131046'; % data timestamp, copy from fig_3_DL_init file

nnidx = [1 4];          % dict sizes to use for the current plot

% Dictionary update routines
updates = {'aksvd', 'nsgk', 'pnsgk'};
% Unused atoms replacement strategy
replatom = 'random';

% Data output
datadir = 'data\';
dataprefix = 'fig_3_11_laols';

% Training images
imdir = 'img\';
images = {'barbara.png', 'boat.png', 'house.png', 'lena.png', 'peppers.png'};
%%-------------------------------------------------------------------------
addpath(genpath('DL'));

timestamp = datestr(now, 'yyyymmddHHMMss');

load([datadir 'fig_3_DL_init-' ots '.mat'], 'Y', 'D0');

methods = length(updates);

% Header
fprintf('Basic Algorithms');
fprintf('\n\tImages: ');
for i = 1:length(images)
    fprintf('%s ', images{i});
end
fprintf('\n\tMethods: ');
for i = 1:methods
    fprintf('%s ', updates{i});
end
fprintf('\n\tParameters: m=%d N=%d K=%d rounds=%d replatoms=%s', ...
    m, N, K, rounds, replatom);

for i = nnidx
    n = nn(i);
for s = ss
    fprintf('\n(n=%d,s=%d): ', n, s); 
    Dall = zeros(rounds,methods, m, n);
    Xall = zeros(rounds,methods, n, N);
    errs = zeros(rounds,methods, K);
    criteria = zeros(rounds,methods, K);
    times = zeros(rounds,methods);    
    for r = 1:rounds
        fprintf('%d', mod(r, 10)); 
        
        %% Rounds
        Yr = zeros(m,N);
        D0r = zeros(m,n);       
        D0r(:,:) = D0{r,i};
        Yr(:,:) = Y{r};

        for j = 1:methods
            fprintf('%s', updates{j}(1));
            time_start = clock;
            [Dall(r,j,:,:), Xall(r,j,:,:), errs(r,j,:)] = ...
                DL(Yr, D0r, s, K, str2func(updates{j}), ...
                'replatom', replatom, 'spfunc', @laols, 'spopts', Lahead);
            time_end = clock;
            times(r,j) = etime(time_end,time_start);  
        end
    end
    %% Write out data
    matfile = sprintf('%s%s-m%d-n%d-N%d-s%d-K%d-%s.mat', ...
         datadir, dataprefix, m, n, N, s, K, timestamp);    
    save(matfile,'updates', 'replatom', 'm', 'n', 'N', 's', 'K', ...
        'Dall', 'Xall', 'errs', 'criteria','times');
end % sparsity loop
end % atoms loop