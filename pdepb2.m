% function [rhs, cu, cux, cuxx, cuy, cuxy, cuyy, cuz, cuxz, cuyz, cuzz] ...
%     = pdepb(x, y, z)
%
% returns the values of the PDE coefficient functions and right side
% ut - Lu = g % z plays the role of t

function [rhs, cu, cux, cuxx, cuy, cuxy, cuyy] = pdepb2(x, y, z)

global Uno Uname;
global etaA etaB etaC nu mu R eee;

oo = ones(size(x)); zz = zeros(size(x));
% watch out!!! -- parabolic
cuz = oo; cuxz = zz; cuyz = zz; cuzz = zz;

if length(y) == 1, y = repmat(y, size(x)); end

%switch Uno
%case {-2}
%otherwise

cu   =  x + y + 1/(z+1);
cux  =  sin(x+2*y) + sin(z);
cuxx =  exp(2*x-y) + cos(z);
cuy  =  cos(2*x+y)*exp(-z);
cuxy = zz; % not implemented yet -- rhscfdpb2
cuyy =  1 + x.*y + 2/(z+1);

%end

[u, ux, uxx, uy, uxy, uyy, uz, uxz, uyz, uzz] = truevd3(x, y, z);

rhs = uz - (cu.*u + cux.*ux + cuxx.*uxx + cuy.*uy + cuxy.*uxy + cuyy.*uyy);
