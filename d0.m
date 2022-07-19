function [r] = d0(v)
% make a sparse diagonal matrix with v the main diagonal
m = length(v(:));
r = spdiags(v(:),0,m,m);
end % d0

