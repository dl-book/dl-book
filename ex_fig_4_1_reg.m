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

%% Generate Figure 4.1: RMSE for standard and regularized DL
close all; fclose all; format compact;
%%-------------------------------------------------------------------------
m = 16;         % problem dimension
n = 32;         % number of atoms in the dictionary
NN = 64:32:512; % number of training signals
s = 4;          % sparsity constraint
reg = 0.1;      % regularization
iters = 50;     % DL iterations
rounds = 50;    % test rounds

% Data location
datadir = 'data\';
dataprefix = 'fig_4_1_reg';
ts = '20170818164513';  % original timestamp

% Dictionary update routines
updates = {'ksvd', 'ksvd_reg', 'aksvd', 'aksvd_reg', 'simco', 'simco_reg'};
color = ['r', 'k', 'm', 'b', 'g', 'y'];

% Figures output
figdir = 'fig\';
figprefix = 'fig_4_1_reg';
%%-------------------------------------------------------------------------
methods = length(updates);

%% Fetch data
allerrs = zeros(length(NN), rounds, methods);
for i = 1:length(NN)
    N = NN(i);
    matfile = sprintf('%s%s-m%d-n%d-N%d-s%d-i%d-%s.mat', ...
         datadir, dataprefix, m, n, N, s, iters, ts);
    load(matfile, 'errs');
    allerrs(i,:,:) = min(errs,[],3);  
end

%% Draw figure
rmsedata = squeeze(mean(allerrs, 2));
f = figure(1);
for up = 1:methods
    plot(NN, squeeze(rmsedata(:,up)), color(up), 'Linewidth', 1);
    xlim([NN(1) NN(end)]);
    hold on;
end
grid on;
hold off;
xticks(NN);
set(gca,'TickLabelInterpreter','latex');
lgd = legend({'K-SVD', 'K-SVDr', 'AK-SVD', 'AK-SVDr', 'SimCO', 'SimCOr'}, ...
    'FontSize', 12, 'Location', 'SouthEast');
xlabel('$N$', 'interpreter', 'latex');
ylabel('RMSE', 'interpreter', 'latex');