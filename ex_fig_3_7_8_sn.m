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

%% Basic Algorithms: DL on images -- generate atom and sparisty figures
%% run_DL.m needs to be executed before this
clear; clc; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
m = 64;                 % problem dimension
nn = 128:64:512;        % number of atoms in the dictionary
N = 4000;               % number of training signals
ss = 4:2:12;            % sparsity constraint
K = 500;                % DL iterations
rounds = 10;            % test rounds
ots = '20170829131046'; % data timestamp, copy from fig_3_DL file

nnidx = [1 4];          % dict sizes to use for the current plot
ssidx = [2 5];          % sparisty levels to use for the current plot

% Dictionary update routines
updates = {'MOD', 'sgk', 'ksvd', 'aksvd', 'nsgk', 'paksvd', 'pnsgk'};
% Unused atoms replacement strategy
replatom = 'random';

% Data output
datadir = 'data\';
dataprefix = 'fig_3_DL';

% Figures output
figdir = 'fig\';   %racheta
figprefix = 'fig_3_78_sn';

color = ['y', 'm', 'c', 'r', 'g', 'b', 'k'];

% Training images
imdir = 'img\'; %racheta
images = {'barbara.png', 'boat.png', 'house.png', 'lena.png', 'peppers.png'};
%%-------------------------------------------------------------------------

methods = length(updates);

subplot(2,1,1)
iplot = 1;
f = figure(1);
%% Gather sparsity data
allerrs = zeros(length(nn), length(ss), rounds, methods, K);
for i = nnidx
    n = nn(i);
    for j = 1:length(ss)
        s = ss(j);
        %% Read data
        matfile = sprintf('%s%s-m%d-n%d-N%d-s%d-K%d-%s.mat', ...
             datadir, dataprefix, m, n, N, s, K, ots);    
        load(matfile,'errs');
        allerrs(i,j,:,:,:) = errs;
    end
end
rmsedata = squeeze(mean(allerrs,3));
nsdata = min(rmsedata,[],4);

for i = nnidx
    n = nn(i);
    figure(1);
    subplot(2,1,iplot);
    for up = 1:methods
        plot(ss, nsdata(i,:,up), color(up), 'Linewidth', 1);
        hold on;
    end
    set(gca, 'Xtick', ss); xlim([ss(1) ss(end)]);
    text(0.45,0.9,['$n=' num2str(n) '$'],...
          'units','normalized','verticalalignment','top', ...
          'Interpreter','LaTex')
    grid on;
    hold off;
    iplot = iplot+1;
end
set(gca,'TickLabelInterpreter','latex');
lgd = legend({'MOD', 'SGK', 'K-SVD', 'AK-SVD', 'NSGK', 'PAK-SVD', 'P-NSGK'});
clear rmsedata nsdata;


subplot(2,1,1)
iplot = 1;
f = figure(2);
%% Gather atom data
allerrs = zeros(length(nn), length(ss), rounds, methods, K);
for i = 1:length(nn)
    n = nn(i);
    for j = ssidx
        s = ss(j);
        %% Read data
        matfile = sprintf('%s%s-m%d-n%d-N%d-s%d-K%d-%s.mat', ...
             datadir, dataprefix, m, n, N, s, K, ots);    
        load(matfile,'errs');
        allerrs(i,j,:,:,:) = errs;
    end
end
rmsedata = squeeze(mean(allerrs,3));
nsdata = min(rmsedata,[],4);

for i = ssidx
    s = ss(i);
    figure(2);
    subplot(2,1,iplot);
    for up = 1:methods
        plot(nn, nsdata(:,i,up), color(up), 'Linewidth', 1);
        hold on;
    end
    set(gca, 'Xtick', nn);xlim([nn(1) nn(end)]);
    text(0.45,0.9,['$s=' num2str(s) '$'],...
          'units','normalized','verticalalignment','top', ...
          'Interpreter','LaTex')
    grid on;
    hold off;
    iplot = iplot+1;
end

set(gca,'TickLabelInterpreter','latex');
lgd = legend({'MOD', 'SGK', 'K-SVD', 'AK-SVD', 'NSGK', 'PAK-SVD', 'P-NSGK'});
clear rmsedata nsdata;