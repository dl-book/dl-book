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

% Example described in Section 9.6.2: 
% Nystrom method with simple two-class problem

% BD 23.12.2017

type_clas = 0;  % classifier used for the actual classification
                % 0 - kernel discriminative DL
                % 1 - SRC with kernel DL

% general parameters
m = 2;
c = 2;         % number of classes
Nc = 100;      % number of training signals per class
N = c*Nc;
n = 60;        % number of atoms
nv = [20; 20]; % number of atoms for SRC
s = 1;         % sparsity level
alpha = 0.1;   % trade-off parameter in discriminative DL objective
nstd = 0.1;    % noise power (for training signals)
nstdt = 0.25;  % noise power for test signals

% kernel parameters
ker_code = 0;
switch ker_code
  case 0,     % RBF
    sigma = 1;
    mykernel = @(x,y) exp( - (x-y)'*(x-y)/2/sigma/sigma);
  case 1,     % polynomial
    a = 1;
    d = 2;
    mykernel = @(x,y) (x'*y + a)^d;
end

% Nystrom parameters
mhat = 30;
p = min(200,N);

% create classes
Y = zeros(m,N);
H = zeros(c,N);
for i = 1:Nc
  Y(:,i) = nstd*randn(m,1);
  v = randn(m,1);
  v = v/norm(v);  % norm 1, random direction
  Y(:,Nc+i) = v*(1 + nstd*randn);
  H(1,i) = 1;
  H(2,Nc+i) = 1;
end

% Nystrom matrix
C = zeros(N,p);
cols = randperm(N,p);
for i = 1:N
  for j = 1:p
    C(i,j) = mykernel(Y(:,i), Y(:,cols(j)));
  end
end
[Q,D] = eig((C(1:p,1:p)+C(1:p,1:p)')/2);
[d,ii] = sort(diag(D), 'descend');
Qhat = Q(:,1:mhat) ./ repmat(sqrt(d(1:mhat))',p,1);
Yhat = Qhat'*C';

% train classifier
switch type_clas
  case 0,       % discriminative DL
    [W, D] = clas_discrim_dl(Yhat, H, n, s, alpha);
  case 1,       % SRC with DL
    D = clas_src_dl(Yhat, H, nv, s);
    n = sum(nv);
end

% test classifier
Ntc = 2000;       % number of test signals for each class
Nt = c*Ntc;
Yt = zeros(m,Nt);
for i = 1:Ntc  % generate signals
  Yt(:,i) = nstdt*randn(m,1);
  v = randn(m,1);
  v = v/norm(v);  % norm 1, random direction
  Yt(:,Ntc+i) = v* (1+nstdt*randn);
end

Ythat = zeros(mhat,Nt);
for i = 1:Nt
  v = zeros(p,1);
  for j = 1:p
    v(j) = mykernel(Y(:,cols(j)), Yt(:,i));
  end
  Ythat(:,i) = Qhat' * v;
end

X = omp(Ythat, D, s);
switch type_clas
  case 0,       % discriminative DL
    Ht = W*X;   % labels for test data
    c1 = find(Ht(1,:) > Ht(2,:));
  case 1,       % SRC with DL
    c1 = find(sum(X(1:nv(1),:)));
end
c2 = setdiff(1:Nt, c1);

clas_succ = (length(find(c1<=Ntc)) + length(find(c2>Ntc)) ) / Nt * 100

% plot training points
figure(1)
plot(Y(1,1:Nc), Y(2,1:Nc), 'or')
hold on
plot(Y(1,Nc+1:N), Y(2,Nc+1:N), 'og')
hold off
grid
%print -depsc fig_nystrom_train.eps

% plot test points
figure(2)
plot(Yt(1,c1), Yt(2,c1), '.r')
hold on
plot(Yt(1,c2), Yt(2,c2), '.g')
hold off
grid
axis([-1.5 1.5 -1.5 1.5])
%print -depsc fig_nystrom_test_1.eps
