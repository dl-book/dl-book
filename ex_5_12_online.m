% Copyright (c) 2017 Bogdan Dumitrescu <bogdan.dumitrescu@acse.pub.ro>
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

% Example 5.12. Compare online DL algorithms (online and RLS DL)
% Generates parts of Figure 5.1 (depending on the choice of the forgetting factor)

% BD 26.05.2017

m = 16;         % problem dimension
n = 40;         % number of atoms in the dictionary
N = 1000;       % number of training signals
s = 4;          % sparsity constraint
ff = 0.999;     % forgetting factor
snr = 20;       % signal-to-noise ratio in dB
iternum = 40*N; % number of iterations
Nruns = 10;

err1 = zeros(iternum,1);  % average RMSE for online DL
err2 = zeros(iternum,1);  % average RMSE for RLS DL

%%-------------------------------------------------------------------------

for k = 1 : Nruns
  k
  % Generate data
  [~, ~, Y] = gen_synth_data(m,n,N,s,snr);
  D0 = normc(randn(m,n));

  % online DL
  [D1, X1, err] = online_dl(Y, D0, s, iternum, ff);
  err1 = err1 + err.*err;

  % RLS DL
  [D2, X2, err] = rls_dl(Y, D0, s, iternum, ff);
  err2 = err2 + err.*err;
end

err1 = sqrt(err1/Nruns);
err2 = sqrt(err2/Nruns);

plot(err1)
hold on
plot(err2,'r')
hold off
grid
csize = 20;
xlabel('iterations', 'interpreter', 'latex', 'FontSize', csize );
ylabel('RMSE', 'interpreter', 'latex', 'FontSize', csize );
h = legend('Online-CD DL', 'RLS-DL');
set(h, 'interpreter', 'latex', 'FontSize', csize );

%print -depsc fig_online.eps

