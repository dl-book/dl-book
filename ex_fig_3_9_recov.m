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

%% Basic Algorithms: DL on synthetic data -- subplot book figure
%% run_recov.m needs to be executed before this
clear; clc; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
m = 64;                 % problem dimension
nn = [128 320];         % number of atoms in the dictionary
N = 4000;               % number of training signals
ss = [6 12];            % sparsity constraint
KK = [150 300];         % max number of iterations to plot
K = 500;                % DL iterations
snrs = Inf;             % Signal to Noise Ratio
rounds = 10;            % test rounds
ots = '20170904120759'; % data timestamp, copy from fig_3_9_recov file

% Dictionary update routines
updates = {'MOD', 'sgk', 'ksvd', 'aksvd', 'nsgk', 'paksvd', 'pnsgk'};

% Data output
datadir = 'data\';
dataprefix = 'fig_3_9_recov';

% Figures output
figdir = 'fig\';
figprefix = 'fig_3_9_recov';

color = ['y', 'm', 'c', 'r', 'g', 'b', 'k'];

%%-------------------------------------------------------------------------

methods = length(updates);
nsnrs = length(snrs);

% Gather error data
allerrs = zeros(length(nn), length(ss), nsnrs, rounds, methods, K);
allrecov = zeros(length(nn), length(ss), nsnrs, rounds, methods, K);
for i = 1:length(nn)
    n = nn(i);
    for j = 1:length(ss)
        s = ss(j);
        for k = 1:nsnrs
            isnr = snrs(k);
            %% Read data
            matfile = sprintf('%s%s-m%d-n%d-N%d-s%d-K%d-snr%d-%s.mat', ...
             datadir, dataprefix, m, n, N, s, K, isnr, ots);
            load(matfile,'errs','recov');
            allerrs(i,j,k,:,:,:) = errs;
            allrecov(i,j,k,:,:,:) = recov;
        end
    end
end

% Build figures
rmsedata = squeeze(mean(allerrs,4));
nsdata = min(rmsedata,[],5);

recovdata = squeeze(mean(allrecov,4));
percentdata = (recovdata./nn')*100;
recovnsdata = recovdata(:,:,:,:,end);
percentnsdata = (recovnsdata./nn')*100;

subplot(2,2,1)
splot = [];
iplot = 1;
f = figure(1);
for i = 1:length(nn)
    n = nn(i);
    for j = 1:length(ss)
        s = ss(j);
        Kmax = KK(j);
        for k = 1:nsnrs
            isnr = snrs(k);

            figure(1);
            splot(iplot) = subplot(2,2,iplot);
            set(gca,'TickLabelInterpreter','latex');
            for up = 1:methods
                plot(1:Kmax, squeeze(percentdata(i,j,up,1:Kmax)), color(up), 'Linewidth', 1);
                hold on;
            end
            text(0.4,0.8,['$n=' num2str(n) '$ $s=' num2str(s) '$'],...
                'units','normalized','verticalalignment','top', ...
                'Interpreter','LaTex')           
            grid on;
            hold off;
            iplot = iplot+1;
        end
    end
end

set(gca,'TickLabelInterpreter','latex');
lgd = legend({'MOD', 'SGK', 'K-SVD', 'AK-SVD', 'NSGK', 'PAK-SVD', 'P-NSGK'});