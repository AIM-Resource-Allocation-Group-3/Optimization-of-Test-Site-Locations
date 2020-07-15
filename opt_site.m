function[x] = opt_site(m,n,currentI)
% f - scoring equation

I = currentI; %number of infected individuals at each node
C = 1; %number of pop-up testing sites
R = zeros(m+2,n+2); %number of fixed testing sites at location (x,y)
% BC
R(1,:) = 0; R(end,:) = 0; R(:,1) = 0; R(:,end) = 0;
R(3,3) = 1;

% fixed testing location modifier
syms M(b)
M(b) = piecewise(b==1,0.2,b==2,0.1,b==3,0.05);
% pop-up testing location modifier
syms G(b)
G(b) = piecewise(b==1,0.2,b==2,0.1,b==3,0);

% Objective function that does not include pop-up modifier
L = zeros(m,n);
for x = 2:m+1
    for y = 2:n+1
        L(x-1,y-1) = I(x-1,y-1)*(1 - M(1)*R(x,y) - M(2)*(R(x+1,y) + R(x-1,y) +...
        R(x,y-1) + R(x,y+1)) - M(3)*(R(x+1,y+1) + R(x-1,y-1) +...
        R(x+1,y-1) + R(x-1,y+1)));
    end
end

intcon = 1:m*n;
A = [];
b = [];
Aeq = ones(1,m*n);
beq = [C];
lb = zeros(m*n,1);
ub = ones(m*n,1);
[x,fval] = intlinprog(-L(:),intcon,A,b,Aeq,beq,lb,ub);
end