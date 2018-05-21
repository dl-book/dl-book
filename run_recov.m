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

%% Basic Algorithms: Dictionary Recovery
clear; clc; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
m = 64;                 % problem dimension
nn = [128 320];         % number of atoms in the dictionary
N = 4000;               % number of training signals
ss = [6 12];            % sparsity constraint
K = 500;                % DL iterations
snrs = Inf;             % Signal to Noise Ratio
rounds = 10;            % test rounds

% Dictionary update routines
updates = {'MOD', 'sgk', 'ksvd', 'aksvd', 'nsgk', 'paksvd', 'pnsgk'};
% Unused atoms replacement strategy
replatom = 'random';

% Data output
datadir = 'data\';   %racheta
dataprefix = 'fig_3_9_recov';

%%-------------------------------------------------------------------------
addpath(genpath('DL'));

timestamp = datestr(now, 'yyyymmddHHMMss');
Y = cell(rounds,length(nn),length(ss),length(snrs));
D0 = cell(rounds,length(nn));
Dtrue = cell(rounds, length(nn));
for r = 1:rounds
    for i = 1:length(nn)
        [D0{r,i},Dtrue{r,i},~] = ...
            gen_synth_data(m,nn(i),N,ss(1),snrs(1));
        for j = 1:length(ss)
            for k = 1:length(snrs)
                Y{r,i,j,k} = ...
                    gen_synth_sigs(Dtrue{r,i},m,nn(i),N,ss(j),snrs(k));
            end
        end
    end
end
save([datadir 'fig_3_9_recov_init-' timestamp '.mat'], 'Y', 'D0', 'Dtrue');

methods = length(updates);
nsnrs = length(snrs);

% Header
fprintf('Basic Algorithms');
fprintf('\n\tRecovery: ');
for i = 1:length(snrs)
    fprintf('%ddB ', snrs(i));
end
fprintf('\n\tMethods: ');
for i = 1:methods
    fprintf('%s ', updates{i});
end
fprintf('\n\tParameters: m=%d N=%d K=%d rounds=%d replatoms=%s', ...
    m, N, K, rounds, replatom);

for i = 1:length(nn)
    n = nn(i);
for j = 1:length(ss)
    s = ss(j);
    %K = 9*s^2;
for k = 1:nsnrs
    isnr = snrs(k);
    fprintf('\n(n=%d,s=%d,snr=%d): ', n, s, isnr);
    Dall = zeros(rounds,methods, m, n);
    Xall = zeros(rounds,methods, n, N);
    errs = zeros(rounds,methods, K);
    recov = zeros(rounds,methods,K);
    times = zeros(rounds,methods);
    for r = 1:rounds
        fprintf('%d', mod(r, 10)); 
        
        %% Rounds
        Yr = zeros(m,N);
        D0r = zeros(m,n);

        D0r(:,:) = D0{r,i};
        Yr(:,:) = Y{r,i,j,k};

        for up = 1:methods
            fprintf('%s', updates{up}(1));
            erropts = {'recov', 'Dtrue', Dtrue{r,i}};
            time_start = clock;
            [Dall(r,up,:,:), Xall(r,up,:,:), errs(r,up,:), recov(r,up,:)] = ...
                DL(Yr, D0r, s, K, str2func(updates{up}), ...
                'replatom', replatom, 'erropts', erropts);
            time_end = clock;
            times(r,up) = etime(time_end,time_start);
        end
    end
    %% Write out data
    matfile = sprintf('%s%s-m%d-n%d-N%d-s%d-K%d-snr%d-%s.mat', ...
         datadir, dataprefix, m, n, N, s, K, isnr, timestamp);    
    save(matfile,'updates', 'replatom', 'm', 'n', 'N', 's', 'K', ...
        'Dall', 'Xall', 'errs', 'recov', 'times');
end % snr loop
end % sparsity loop
end % atoms loop