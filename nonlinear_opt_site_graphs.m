% Very slow and puts testing sites on boundaries
function[x,fval,exitFlag,Output] = nonlinear_opt_site_graphs(m,n,currentI,neib,dneib)
I = currentI; %number of infected individuals at each node

C = 1; %number of pop-up testing sites
R = zeros(m*n,1); %number of fixed testing sites at each node
% fixed testing location modifier
syms M(b)
M(b) = piecewise(b==1,0.2,b==2,0.1,b==3,0.05);
% pop-up testing location modifier
syms G(b)
G(b) = piecewise(b==1,0.2,b==2,0.1,b==3,0);

M = [M(1), M(2), M(3)];
G = [G(1), G(2), G(3)];

Aeq = [];
beq = [];
A = ones(1,m*n);
b = C;
lb = zeros(m*n,1);
ub = ones(m*n,1);
nonlcon = [];
intcon = 1:m*n;

fun = @(F) objfun(m,n,I,M,R,G,F,neib,dneib);
% options for ga() solver (go to Input->options for full list): https://www.mathworks.com/help/gads/ga.html
    % PlotFcn - displays plot as the optimization is working
    % MaxStallGenerations - one possible stopping criteria (looks at avg relative change)
    % MaxGerations - number of iterations before ga stops
options = optimoptions('ga','PopulationSize',500,'PlotFcn',{@gaplotbestf,@gaplotstopping},...
    'MaxStallGenerations', 300,'MaxGenerations',200);
options = optimoptions(options,'UseVectorized',true);
[x,fval,exitFlag,Output] = ga(fun,m*n,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options);
end

function L = objfun(m,n,I,M,R,G,F,neib,dneib)
    L = zeros(m*n,1);
    for node = 1:m*n
        nnode = neib(node,:); nnode = nnode(nnode~=0);
        dnode = dneib(node,:); dnode = dnode(dnode~=0);
        L(node) = I(node)*(1 - M(1)*R(node) - G(1)*F(node) - M(2)*sum(R(nnode)) ...
            - G(2)*sum(F(nnode)) - M(3)*sum(R(dnode)) - G(3)*sum(F(dnode)));
    end
    L = -sum(F*L);
end