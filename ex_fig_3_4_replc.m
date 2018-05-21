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


%% Basic Algorithms: DL atom replacement figure
%% run_replc.m needs to be executed before this
clear; clc; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
m = 64;                 % problem dimension
nn = [128 192 256 512]; % number of atoms in the dictionary
N = 2048;               % number of training signals
s = 8;                  % sparsity constraint
K = 50;                 % DL iterations
rounds = 10;            % test rounds
ots = '20171120205154'; % data timestamp, copy from fig_3_4 file

% Dictionary update routines
updates = {'aksvd'};
% Unused atoms replacement strategy
repl = {'no', 'random', 'worst', 'zero'};
extra = {'', '', '', 'worstrep'};

% Data output
datadir = 'data\';   %racheta
dataprefix = 'fig_3_4_replc';

% Figures output
figdir = 'fig\';   %racheta
figprefix = 'fig_3_4_replc';

color = ['k', 'r', 'g', 'b'];

% Training images
imdir = 'img\';
images = {'barbara.png', 'boat.png', 'house.png', 'lena.png', 'peppers.png'};
%%-------------------------------------------------------------------------

methods = length(updates);

%% Gather repl data
allerrs = zeros(length(nn), length(repl), rounds, K);
for i = 1:length(nn)
    n = nn(i);
    for j = 1:length(repl)
        replatom = repl{j};
        %% Read data
        matfile = sprintf('%s%s-m%d-n%d-N%d-s%d-K%d-%s-%s.mat', ...
             datadir, dataprefix, m, n, N, s, K, replatom, ots);    
        load(matfile,'Y','errs');
        Ynorm = cellfun(@(x) norm(x,'fro'), Y);
        allerrs(i,j,:,:) = (squeeze(errs)*sqrt(m*N))./repmat(Ynorm,1,K);
    end
end
rmsedata = squeeze(mean(allerrs,3));

subplot(2,2,1)
iplot = 1;
f = figure(1);
for i = 1:length(nn)
    n = nn(i);
    figure(1);
    subplot(2,2,iplot);
    for j = 1:length(repl)
        plot(1:K, squeeze(rmsedata(i,j,:)), color(j), 'Linewidth', 1);
        hold on;
    end
    %set(gca, 'Xtick', ss); xlim([ss(1) ss(end)]);
    xlim([1 K]);
    text(0.45,0.9,['$n=' num2str(n) '$'],...
          'units','normalized','verticalalignment','top', ...
          'Interpreter','LaTex')
    grid on;
    hold off;
    iplot = iplot+1;
end
set(gca,'TickLabelInterpreter','latex');
lgd = legend({'None', 'Random', 'Worst', 'Post'});