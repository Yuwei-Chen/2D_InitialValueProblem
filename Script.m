%Test function set-up
format compact
global Uno Uname;
Uno = 17; Uname = '';

%Initial Setup
ax = 0; bx = 1; % domain on x
ay = 0; by = 1; % domain on y

ntimes = 5; errg = zeros(1, ntimes); %number of grid sizes to test

T = 1; %time
theta = 1/2; % 0 for forward, 1 for backward, 1/2 for Crank-Nicolson
Nt = 100; ht = T/Nt; % uniform timestep

for nn = 1:ntimes
    n = 2^(nn+2); ngrid = n+1; nint(nn) = n;
    neq = (n-1)^2; numeq = neq;	% (n-1)*(m-1) for D
    
    %set up of non-uniform grid
    hx_u = (bx-ax)/n;
    gridx_u = ax + hx_u*[0:n];
    gridx = gridx_u.^2;
    hx = gridx(2:end) - gridx(1:end-1);
    
    hy_u = (by-ay)/n;
    gridy_u = ax + hx_u*[0:n];
    gridy = gridy_u.^2;
    hy = gridy(2:end) - gridy(1:end-1);
    
    %initial value of u
    u = [];
    for i = 2:(length(gridx)-1)
        for j = 2:(length(gridy)-1)
        [temp] = truevd3(gridx(i), gridy(j), 0);
        u = [u; temp];
        end
    end
    t = 0;
    
    %solving IVP
    [bcv, rhsf, coefs] = rhscfd2d( gridx, gridy, 0 );
    A = cfdmat2d(gridx, gridy, coefs);
    
    for j = 1:Nt
        
        %timestep
        ht = ht;  %temporary little time-step
        t = t + ht;
        
        %backward matrix fot this timestep
        [bcv2, rhsf2, coefs2] = rhscfd2d( gridx, gridy, t );
        A2 = cfdmat2d(gridx, gridy, coefs2);
        
        %setup of system and solve
        rhs = ( speye(size(A))/ht+(1-theta)*A ) * u ...
         + rhsf*(1-theta) + rhsf2*theta + theta*bcv2 + (1-theta)*bcv;
        AA = speye(size(A))/ht - theta*A2;
        u = AA\rhs;
        
        bcv = bcv2; rhsf = rhsf2; coefs = coefs2; A = A2;
    end
    
    %Inf-norm abs error (approximation and true value)
    errg = errorfd(ngrid, gridx, gridy, n, u, T, nn, errg);
    
end

%OUTPUT
[udummy] = truevd3(ax, ay, T);
disp(['U = ' Uname ' = {' num2str(Uno) '}'])
disp(['domain [' num2str(ax)  ', ' num2str(bx) '] X [' num2str(ay) ', '...
    num2str(by) '] at time ' num2str(T)])

nint
format short e
disp('error on grid points')
errg
format short
disp('order of convergence')
errg(:, :) = max(errg(:, :),  0.222044604925e-15);
LogNintRatio = log(nint(1, 2:ntimes)./nint(1, 1:ntimes-1));
LogNintRatioMat = repmat(LogNintRatio, size(errg, 1), 1);
if ntimes > 1
    convg = log(errg(:, 1:ntimes-1)./errg(:, 2:ntimes))./LogNintRatioMat
end