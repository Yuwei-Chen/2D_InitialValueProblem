% function [u, ux, uxx, uy, uxy, uyy, uz, uxz, uyz, uzz] = U(x, y, z)
%
% returns the values of the solution function,
% and derivatives up to 2nd order on x, y, i.e. u, ux, uxx, uy, uxy, uyy,
% uz, uxz, uyz, uzz
% 1: u
% 2: ux
% 3: uxx
% 4: uy
% 5: uxy
% 6: uyy
% 7: uz
% 8: uxz
% 9: uyz
% 10: uzz
% x can be a vector/matrix. If x is a vector/matrix,
% then y must be either scalar or vector/matrix of same size as x.
% z must be scalar

function [u, ux, uxx, uy, uxy, uyy, uz, uxz, uyz, uzz] = truevd3(x, y, z)

global Uno Uname;

oo = ones(size(x)); zz = zeros(size(x));
yy = y;
if length(y) == 1, y = repmat(y, size(x)); end
switch Uno
case {1113}
    Uname ='exp(x).*abs(y-1/2) + z.^2';
    u   = exp(x).*abs(y-1/2) + z.^2;
    ux  = exp(x).*abs(y-1/2);
    uxx = exp(x).*abs(y-1/2);
    uy  = (y>=1/2).*exp(x) - (y<1/2).*exp(x);
    uxy = uy;
    uz  = 0;
    uzz = 2*z;
    uxz = 0;
    uyz = 0;
    uzz = 2;
case {165}
    Uname = 'x.^(13/2) .* y.^(13/2).* z.^(13/2)';
    ux0 = x.^(13/2);      ux1 = 13/2 .* x.^(11/2); ux2 = 143/4 .* x.^(9/2);
    uy0 = y.^(13/2);      uy1 = 13/2 .* y.^(11/2); uy2 = 143/4 .* y.^(9/2);
    uz0 = z.^(13/2);      uz1 = 13/2 .* z.^(11/2); uz2 = 143/4 .* z.^(9/2);
    u  = ux0.*uy0.*uz0;   ux  = ux1.*uy0.*uz0;   uxx = ux2.*uy0.*uz0;
    uy = ux0.*uy1.*uz0;   uxy = ux1.*uy1.*uz0;   uyy = ux0.*uy2.*uz0;
    uz = ux0.*uy0.*uz1;   uxz = ux1.*uy0.*uz1;   uyz = ux0.*uy1.*uz1;
    uzz= ux0.*uy0.*uz2;
case {155}
    Uname = 'x.^(11/2) .* y.^(11/2).* z.^(11/2)';
    ux0 = x.^(11/2);      ux1 = 11/2 .* x.^(9/2); ux2 = 99/4 .* x.^(7/2);
    uy0 = y.^(11/2);      uy1 = 11/2 .* y.^(9/2); uy2 = 99/4 .* y.^(7/2);
    uz0 = z.^(11/2);      uz1 = 11/2 .* z.^(9/2); uz2 = 99/4 .* z.^(7/2);
    u  = ux0.*uy0.*uz0;   ux  = ux1.*uy0.*uz0;   uxx = ux2.*uy0.*uz0;
    uy = ux0.*uy1.*uz0;   uxy = ux1.*uy1.*uz0;   uyy = ux0.*uy2.*uz0;
    uz = ux0.*uy0.*uz1;   uxz = ux1.*uy0.*uz1;   uyz = ux0.*uy1.*uz1;
    uzz= ux0.*uy0.*uz2;
case {145}
    Uname = 'x.^(9/2) .* y.^(9/2).* z.^(9/2)';
    ux0 = x.^(9/2);       ux1 = 9/2 .* x.^(7/2); ux2 = 63/4 .* x.^(5/2);
    uy0 = y.^(9/2);       uy1 = 9/2 .* y.^(7/2); uy2 = 63/4 .* y.^(5/2);
    uz0 = z.^(9/2);       uz1 = 9/2 .* z.^(7/2); uz2 = 63/4 .* z.^(5/2);
    u  = ux0.*uy0.*uz0;   ux  = ux1.*uy0.*uz0;   uxx = ux2.*uy0.*uz0;
    uy = ux0.*uy1.*uz0;   uxy = ux1.*uy1.*uz0;   uyy = ux0.*uy2.*uz0;
    uz = ux0.*uy0.*uz1;   uxz = ux1.*uy0.*uz1;   uyz = ux0.*uy1.*uz1;
    uzz= ux0.*uy0.*uz2;
case {135}
    Uname = 'x.^(7/2) .* y.^(7/2).* z.^(7/2)';
    ux0 = x.^(7/2);      ux1 = 7/2 .* x.^(5/2); ux2 = 35/4 .* x.^(3/2);
    uy0 = y.^(7/2);      uy1 = 7/2 .* y.^(5/2); uy2 = 35/4 .* y.^(3/2);
    uz0 = z.^(7/2);      uz1 = 7/2 .* z.^(5/2); uz2 = 35/4 .* z.^(3/2);
    u  = ux0.*uy0.*uz0;   ux  = ux1.*uy0.*uz0;   uxx = ux2.*uy0.*uz0;
    uy = ux0.*uy1.*uz0;   uxy = ux1.*uy1.*uz0;   uyy = ux0.*uy2.*uz0;
    uz = ux0.*uy0.*uz1;   uxz = ux1.*uy0.*uz1;   uyz = ux0.*uy1.*uz1;
    uzz= ux0.*uy0.*uz2;
case {119}
    Uname ='exp(x+y+z)';
    ux0 = exp(x);         ux1 = exp(x);          ux2 = exp(x);
    uy0 = exp(y);         uy1 = exp(y);          uy2 = exp(y);
    uz0 = exp(z);         uz1 = exp(z);          uz2 = exp(z);
    u  = ux0.*uy0.*uz0;   ux  = ux1.*uy0.*uz0;   uxx = ux2.*uy0.*uz0;
    uy = ux0.*uy1.*uz0;   uxy = ux1.*uy1.*uz0;   uyy = ux0.*uy2.*uz0;
    uz = ux0.*uy0.*uz1;   uxz = ux1.*uy0.*uz1;   uyz = ux0.*uy1.*uz1;
    uzz= ux0.*uy0.*uz2;
    D4X =exp(x+y+z);
case {111}
    Uname ='sin(x).*sin(3*y).*exp(-z)';
    ux0 = sin(x);         ux1 = cos(x);          ux2 =-sin(x);
    uy0 = sin(3*y);       uy1 = 3*cos(3*y);      uy2 =-9*sin(3*y);
    uz0 = exp(-z);        uz1 =-exp(-z);         uz2 = exp(-z);
    u  = ux0.*uy0.*uz0;   ux  = ux1.*uy0.*uz0;   uxx = ux2.*uy0.*uz0;
    uy = ux0.*uy1.*uz0;   uxy = ux1.*uy1.*uz0;   uyy = ux0.*uy2.*uz0;
    uz = ux0.*uy0.*uz1;   uxz = ux1.*uy0.*uz1;   uyz = ux0.*uy1.*uz1;
    uzz= ux0.*uy0.*uz2;
case {110}
    Uname ='sin(x).*sin(y).*exp(-z)';
    ux0 = sin(x);         ux1 = cos(x);          ux2 =-sin(x);
    uy0 = sin(y);         uy1 = cos(y);          uy2 =-sin(y);
    uz0 = exp(-z);        uz1 =-exp(-z);         uz2 = exp(-z);
    u  = ux0.*uy0.*uz0;   ux  = ux1.*uy0.*uz0;   uxx = ux2.*uy0.*uz0;
    uy = ux0.*uy1.*uz0;   uxy = ux1.*uy1.*uz0;   uyy = ux0.*uy2.*uz0;
    uz = ux0.*uy0.*uz1;   uxz = ux1.*uy0.*uz1;   uyz = ux0.*uy1.*uz1;
    uzz= ux0.*uy0.*uz2;
case {39}
    Uname ='x.^4.*y.^4.*z^4';
    ux0 = x.^4;           ux1 = 4*x.^3;          ux2 = 12*x.^2;
    uy0 = y.^4;           uy1 = 4*y.^3;          uy2 = 12*y.^2;
    uz0 = z.^4;           uz1 = 4*z.^3;          uz2 = 12*z.^2;
    u  = ux0.*uy0.*uz0;   ux  = ux1.*uy0.*uz0;   uxx = ux2.*uy0.*uz0;
    uy = ux0.*uy1.*uz0;   uxy = ux1.*uy1.*uz0;   uyy = ux0.*uy2.*uz0;
    uz = ux0.*uy0.*uz1;   uxz = ux1.*uy0.*uz1;   uyz = ux0.*uy1.*uz1;
    uzz= ux0.*uy0.*uz2;
case {29}
    Uname ='x.^3.*y.^3.*z^3';
    ux0 = x.^3;           ux1 = 3*x.^2;          ux2 = 6*x;
    uy0 = y.^3;           uy1 = 3*y.^2;          uy2 = 6*y;
    uz0 = z.^3;           uz1 = 3*z.^2;          uz2 = 6*z;
    u  = ux0.*uy0.*uz0;   ux  = ux1.*uy0.*uz0;   uxx = ux2.*uy0.*uz0;
    uy = ux0.*uy1.*uz0;   uxy = ux1.*uy1.*uz0;   uyy = ux0.*uy2.*uz0;
    uz = ux0.*uy0.*uz1;   uxz = ux1.*uy0.*uz1;   uyz = ux0.*uy1.*uz1;
    uzz= ux0.*uy0.*uz2;
case {28}
    Uname ='x.^3.*y.^3.*z^2';
    ux0 = x.^3;           ux1 = 3*x.^2;          ux2 = 6*x;
    uy0 = y.^3;           uy1 = 3*y.^2;          uy2 = 6*y;
    uz0 = z.^2;           uz1 = 2*z;             uz2 = 2;
    u  = ux0.*uy0.*uz0;   ux  = ux1.*uy0.*uz0;   uxx = ux2.*uy0.*uz0;
    uy = ux0.*uy1.*uz0;   uxy = ux1.*uy1.*uz0;   uyy = ux0.*uy2.*uz0;
    uz = ux0.*uy0.*uz1;   uxz = ux1.*uy0.*uz1;   uyz = ux0.*uy1.*uz1;
    uzz= ux0.*uy0.*uz2;
case {19}
    Uname ='x.x.*y.*y *z *z';
    u  = x.*x.*y.*y *z *z; ux  = 2*x.*y.*y *z *z; uxx = 2    *y.*y *z *z;
    uy = x.*x.*2.*y *z *z; uxy = 2*x.*2.*y *z *z; uyy = x.*x.*2    *z *z;
    uz = x.*x.*y.*y *2 *z; uxz = 2*x.*y.*y *2 *z; uyz = x.*x.*2.*y *2 *z;
    uzz= x.*x.*y.*y *2;
case {18}
    Uname ='(x.*x - x).*(y.*y - y).*(z.*z - z)';
    ux0 = (x.*x - x);     ux1 = 2*x - oo;        ux2 = 2*oo;
    uy0 = (y.*y - y);     uy1 = 2*y - oo;        uy2 = 2*oo;
    uz0 = (z.*z - z);     uz1 = 2*z - oo;        uz2 = 2*oo;
    u  = ux0.*uy0.*uz0;   ux  = ux1.*uy0.*uz0;   uxx = ux2.*uy0.*uz0;
    uy = ux0.*uy1.*uz0;   uxy = ux1.*uy1.*uz0;   uyy = ux0.*uy2.*uz0;
    uz = ux0.*uy0.*uz1;   uxz = ux1.*uy0.*uz1;   uyz = ux0.*uy1.*uz1;
    uzz= ux0.*uy0.*uz2;
case {17}
    Uname ='x.x.*y.*y *z';
    u  = x.*x.*y.*y *z;    ux  = 2*x.*y.*y *z;    uxx = 2    *y.*y *z;
    uy = x.*x.*2.*y *z;    uxy = 2*x.*2.*y *z;    uyy = x.*x.*2    *z;
    uz = x.*x.*y.*y;       uxz = 2*x.*y.*y;       uyz = x.*x.*2.*y;
    uzz= zz;
case {9}
    Uname ='x.*y*z';
    u  = x.*y*z;          ux  = y*z;             uxx = zz;
    uy = x *z;            uxy = z*oo;            uyy = zz;
    uz = x.*y;            uxz = y;               uyz = x;
    uzz= zz;
case {4}
    Uname ='z';
    u  = z;               ux  = zz;              uxx = zz;
    uy = zz;              uxy = zz;              uyy = zz;
    uz = oo;              uxz = zz;              uyz = zz;
    uzz= zz;
case {3}
    Uname ='y';
    u  = y;               ux  = zz;              uxx = zz;
    uy = oo;              uxy = zz;              uyy = zz;
    uz = zz;              uxz = zz;              uyz = zz;
    uzz= zz;
case {1}
    Uname ='x';
    u  = x;               ux  = oo;              uxx = zz;
    uy = zz;              uxy = zz;              uyy = zz;
    uz = zz;              uxz = zz;              uyz = zz;
    uzz= zz;
case {0}
    Uname ='1';
    u  = oo;              ux  = zz;              uxx = zz;
    uy = zz;              uxy = zz;              uyy = zz;
    uz = zz;              uxz = zz;              uyz = zz;
    uzz= zz;
otherwise
    error(['truevd3: no such function ' num2str(Uno)])
end
