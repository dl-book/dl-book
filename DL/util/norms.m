% Copyright (c) 2016 Paul Irofti <paul@irofti.net>
% 
% Public domain.
% 

function n = norms(A)
    n = sqrt(sum(A.^2,1));
end