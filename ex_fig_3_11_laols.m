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


%% Basic Algorithms: DL on images -- generate LAOLS figures
%% run_laols.m and run_DL.m needs to be executed before this
clear; clc; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
m = 64;                 % problem dimension
nn = 128:64:512;        % number of atoms in the dictionary
N = 4000;               % number of training signals
ss = 4:2:12;            % sparsity constraint
Komp = 500;             % DL OMP iterations
Klaols = 100;           % DL LAOLS iterations
rounds = 10;            % test rounds

omptime = '20170829131046';   % omp timestamp (from run_DL)
laolstime = '20170925144730'; % laols timestamp (from run_laols)

% Dictionary update routines
updates = {'MOD', 'sgk', 'ksvd', 'aksvd', 'nsgk', 'paksvd', 'pnsgk'};
laolsupdates = {'aksvd', 'nsgk', 'pnsgk'};
% Unused atoms replacement strategy
replatom = 'random';

% Interesting data
ompmethods = length(updates);
laolsmethods = length(laolsupdates);
nidx = [1 4];           % interesting dictionary sizes
sidx = 2;               % interesting sparsities
upidx = [4 5 7];        % interesting methods

% Data output
datadir = 'data\';
laolsprefix = 'fig_3_11_laols';
ompprefix = 'fig_3_DL';

% Figures output
figdir = 'fig\';
figprefix = 'fig_3_11_laols';

color = ['y', 'm', 'c', 'r', 'g', 'b', 'k'];

% Training images
imdir = 'img\';
images = {'barbara.png', 'boat.png', 'house.png', 'lena.png', 'peppers.png'};
%%-------------------------------------------------------------------------

%% Gather iterations data
% LAOLS
allerrs = zeros(length(nn), length(ss), rounds, laolsmethods, Klaols);
for i = nidx
    n = nn(i);
    for j = sidx
        s = ss(j);
        %% Read data
        matfile = sprintf('%s%s-m%d-n%d-N%d-s%d-K%d-%s.mat', ...
             datadir, laolsprefix, m, n, N, s, Klaols, laolstime);
        load(matfile,'errs');
        allerrs(i,j,:,:,:) = errs;
    end
end
%allerrs = allerrs(nidx,sidx,:,:,:);
laolsdata = squeeze(mean(allerrs,3));
% OMP
allerrs = zeros(length(nn), length(ss), rounds, ompmethods, Komp);
for i = nidx
    n = nn(i);
    for j = sidx
        s = ss(j);
        %% Read data
        matfile = sprintf('%s%s-m%d-n%d-N%d-s%d-K%d-%s.mat', ...
             datadir, ompprefix, m, n, N, s, Komp, omptime);
        load(matfile,'errs');
        allerrs(i,j,:,:,:) = errs;
    end
end
%allerrs = allerrs(nidx,sidx,:,:,:);
ompdata = squeeze(mean(allerrs,3));

K = min(Komp,Klaols);
subplot(2,1,1)
splot = [];
iplot = 1;
f = figure(1);
for i = nidx
    n = nn(i);
    splot(iplot) = subplot(2,1,iplot);
    for j = sidx
        s = ss(j);
        %% First get ylim
        figure(1);
        for k = 1:laolsmethods
            plot(5:K, squeeze(laolsdata(i,j,k,5:K)), ...
                color(upidx(k)), 'Linewidth', 1);
            hold on;
        end
        for k = upidx
            plot(5:K, squeeze(ompdata(i,j,k,5:K)), color(k), 'Linewidth', 1);
            hold on;
        end
        yl = ylim;
        hold off;

        %% Do the real figure
        figure(1);
        plot([0 0],[0 0],'w-'); hold on;% dummy for OMP text
        for k = 1:laolsmethods
            plot(1:K, squeeze(laolsdata(i,j,k,1:K)), ...
                ['--' color(upidx(k))], 'Linewidth', 1);
            hold on;
        end
        plot([0 0],[0 0],'w-'); hold on;% dummy for LAOLS text
        for k = upidx
            plot(1:K, squeeze(ompdata(i,j,k,1:K)), ...
                color(k), 'Linewidth', 1);
            hold on;
        end

        text(0.3,0.9,['$n=' num2str(n) '$'],...
          'units','normalized','verticalalignment','top', ...
          'Interpreter','LaTex')
        ylim(yl);
        set(gca, 'FontSize', 12);
        grid on;
        hold off;
    end
    iplot = iplot+1;
end
lgd = legend({'\textbf{LAOLS}', 'AK-SVD', 'NSGK', 'P-NSGK', ...
    '\textbf{OMP}', 'AK-SVD', 'NSGK', 'P-NSGK'},...
    'FontSize', 12);
set(lgd,'interpreter','latex');