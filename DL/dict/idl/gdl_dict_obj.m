function [f, grad] = gdl_dict_obj(d_vec, X, XXt, YXt, Y, p, n, m, gamma)
% Objective function to minimize in dictionary update of generative
% dictionary learning.
%
% See writing/dictionary_update_with_coherence_penalty.lyx for
% documentation.

D = reshape(d_vec,p,n);

%f = 1/p/m*norm(D*X-Y,'fro')^2+gamma/n^2*norm(D'*D-eye(n),'fro')^2;
f = norm(D*X-Y,'fro')^2+gamma*norm(D'*D-eye(n),'fro')^2;

if nargout > 1
    %Grad = 1/p/m*(2*D*XXt - 2*YXt) + 4*gamma/n^2*((D*D')*D - D);
    Grad = 2*(D*XXt - YXt) + 4*gamma*((D*D')*D - D);
    grad = Grad(:);
end
end
