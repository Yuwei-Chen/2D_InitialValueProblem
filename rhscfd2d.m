function [bcv, rhsf, coefs ] = rhscfd2d( gridx, gridy, t )
%RHSCFD
%   Dirichlet conditions

n = length(gridx) - 1; % number of grid on x
m = length(gridy) - 1; % number of grid on y
neq = (n - 1)*(m - 1);	% assume Dirichlet conditions
hx = gridx(2:end) - gridx(1:end-1);
hy = gridy(2:end) - gridy(1:end-1);

%boundary condition and modification vector
for i = 1:(n-1)
    bv1(i) = bc2(gridx(end), gridy(i+1), t,1);
    bv3(i) = bc2(gridx(1), gridy(i+1), t, 3);
end
for i = 1:(m-1)
    bv2(i) = bc2(gridx(i+1), gridy(1),t, 2);
    bv4(i) = bc2(gridx(i+1), gridy(end),t, 4);
end
bv5 = bc2(gridx(end), gridy(1),t,5); bv6 = bc2(gridx(1), gridy(1),t,6);
bv7 = bc2(gridx(1), gridy(end),t,7); bv8 = bc2(gridx(end), gridy(end),t,8);
bcv = zeros(neq, 1); 

%%RHS Source Functions and COEFs
counter = 1;
for i = 1:(n-1)
    for j = 1:(m-1)
        px = gridx(i+1); py = gridy(j+1);
        [rhsf(counter, 1), coefs(counter, 1), coefs(counter, 2), ...
            coefs(counter, 3), coefs(counter, 4), coefs(counter, 5),...
            coefs(counter, 6)] = pdepb2(px, py, t);
        
        %u_xx and u_x
        if i == 1
            bcv(counter) = bcv(counter) + bv3(j)*...
                (2*coefs(j,3) - coefs(j,2)*hx(2))/(hx(1)*(hx(1)+hx(2)));
        elseif i == n-1
            bcv(counter) = bcv(counter) + bv1(j)*...
                (2*coefs(counter,3) + coefs(counter,2)*hx(end-1))/(hx(end)*(hx(end)+hx(end-1)));
        end
        
        %u_yy and u_y
        if j == 1
            bcv(counter) = bcv(counter) + bv2(i)*...
                (2*coefs(counter,6) - coefs(counter,4)*hy(2))/(hy(1)*(hy(1)+hy(2)));
        elseif j == m-1
            bcv(counter) = bcv(counter) + bv4(i)*...
                (2*coefs(counter,6) + coefs(counter,4)*hy(end-1))/(hy(end)*(hy(end)+hy(end-1)));
        end
        
        %u_xy
        if i == 1 && j == 1
            bcv(counter) = bcv(counter) + coefs(counter, 5)*...
                (bv6-bv3(2)-bv2(2))/(hx(1)+hx(2))/(hy(1)+hy(2));
        elseif i == 1 && j == (m-1)
            bcv(counter) = bcv(counter) +  coefs(counter, 5)*...
                (bv3(end-1)+bv4(2)-bv7)/(hx(1)+hx(2))/(hy(end)+hy(end-1));
        elseif i == (n-1) && j == 1
            bcv(counter) = bcv(counter) + coefs(counter, 5)*...
                (bv2(end-1)+bv1(2)-bv5)/(hy(1)+hy(2))/(hx(end)+hx(end-1));
        elseif i == (n-1) && j == (m-1)
            bcv(counter) = bcv(counter) + coefs(counter, 5)*...
                (bv8-bv4(end-1)-bv1(end-1))/(hx(end)+hx(end-1))/(hy(end)+hy(end-1));
        elseif i == 1 && j ~= 1 && j ~= (m-1)
            bcv(counter) = bcv(counter) + coefs(counter, 5)*...
                (bv3(j-1)-bv3(j+1))/(hx(1)+hx(2))/(hy(j)+hy(j+1));
        elseif i == (n-1) && j ~= 1 && j ~= (m-1)
            bcv(counter) = bcv(counter) + coefs(counter, 5)*...
                (-bv1(j-1)+bv1(j+1))/(hx(end)+hx(end-1))/(hy(j)+hy(j+1));
        elseif j == 1 && i ~= 1 && i ~= (n-1)
            bcv(counter) = bcv(counter) + coefs(counter, 5)*...
                (bv2(i-1)-bv2(i+1))/(hx(i)+hx(i+1))/(hy(1)+hy(2));
        elseif j == (m-1) && i ~= 1 && i ~= (n-1)
            bcv(counter) = bcv(counter) + coefs(counter, 5)*...
                (-bv4(i-1)+bv4(i+1))/(hx(i)+hx(i+1))/(hy(end-1)+hy(end));
        end
        
        counter = counter + 1;
    end
end

end