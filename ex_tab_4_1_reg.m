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

%% Generate Table 4.1: RMSE for standard and regularized DL
clear; clc; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
m = 64;                         % problem dimension
n = 256;                        % number of atoms in the dictionary
N = 4096;                       % number of training signals
ss = [4 6 8 10 12];             % sparsity constraint
reg = 0.05;                     % regularization
vanish = 1;                     % regularization vanishing factor
pa = n;                         % number of atoms to update in parallel
iters = 50;                     % DL iterations
rounds = 10;                    % test rounds
itN = 1;                        % SimCO internal iterations
iitN = [10, 20, 30, 40, 50];    % SimCO internal iterations

datadir = 'data\';
dataprefix = 'tab_4_1_reg';
ts = '20170818164513'; % original timestamp

% Dictionary update routines
updates = {'KSVD', 'KSVDr', 'AK-SVD', 'AK-SVDr', 'SimCO(1)', 'SimCOr(1)'};

%%-------------------------------------------------------------------------
methods = length(updates);
minerr = ones(1, length(ss));

fprintf('s\t\t\t%s\n', sprintf('%d\t\t\t', ss));
for j = 1:methods
    fprintf('%s\t', updates{j});
    for s = ss
        matfile = sprintf('%s%s-m%d-n%d-N%d-s%d-i%d-%s.mat', ...
                datadir, dataprefix, m, n, N, s, iters, ts);
        load(matfile, 'errs');
        avg = mean(min(errs(:,j,:),[],3));
        fprintf('\t%f', avg);
        if minerr(find(s == ss)) > avg
            minerr(find(s == ss)) = avg;
        end
        if s ~= ss(end)
            fprintf(' ');
        end
        clear('errs');
    end
    fprintf('\n');
end

for j = 1:length(iitN)
    fprintf('SimCOr(%d)\t', iitN(j));
    for s = ss
        matfile = sprintf('%s%s-itN%d-m%d-n%d-N%d-s%d-i%d-%s.mat', ...
                datadir, dataprefix, iitN(j),m, n, N, s, iters, ts);
        load(matfile, 'errs');
        avg = mean(min(errs(:,1,:),[],3));
        fprintf('\t%f', avg);
        if minerr(find(s == ss)) > avg
            minerr(find(s == ss)) = avg;
        end
        if s ~= ss(end)
            fprintf(' ');
        end
        clear('errs');
    end
    fprintf('\n');
end
fprintf('\n');