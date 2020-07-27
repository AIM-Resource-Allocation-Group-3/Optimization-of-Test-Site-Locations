% Very slow and puts testing sites on boundaries
function[x,fval,exitFlag,Output,pop,scores] = nonlinear_opt_site_graphs(currentI,didt,neib,dneib,possTest,dist,fixneib,fixdneib,possneib,possdneib)
I = currentI; %number of infected individuals at each node
C = 5; %number of pop-up testing sites
numPossTest = length(possTest);
% R = zeros(length(I),1); %number of fixed testing sites at each node
% fixed testing location modifier
syms M(b)
M(b) = piecewise(b==1,3,b==2,2,b==3,1);
% pop-up testing location modifier
syms G(b)
G(b) = piecewise(b==1,3,b==2,2,b==3,1);

M = [M(1), M(2), M(3)];
G = [G(1), G(2), G(3)];

Aeq = [];
beq = [];
A = ones(1,numPossTest);
b = C;
lb = zeros(numPossTest,1);
ub = ones(numPossTest,1);
% ub(possTest) = 1;
nonlcon = [];
intcon = 1:numPossTest;
initPop = zeros(100,numPossTest);

fun = @(F) objfun(I,M,G,F,didt,neib,dneib,possTest,dist,fixneib,fixdneib,possneib,possdneib);
% options for ga() solver (go to Input->options for full list): https://www.mathworks.com/help/gads/ga.html
    % PlotFcn - displays plot as the optimization is working
    % MaxStallGenerations - one possible stopping criteria (looks at avg relative change)
    % MaxGerations - number of iterations before ga stops
options = optimoptions('ga','PopulationSize',400,'PlotFcn',{@gaplotbestf,@gaplotstopping},...
    'FunctionTolerance',10e-3,'MaxStallGenerations',15,'MaxGenerations',50,'InitialPopulationMatrix',initPop);
options = optimoptions(options,'UseVectorized',true);
[x,fval,exitFlag,Output,pop,scores] = ga(fun,numPossTest,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options);
end

function L = objfun(I,M,G,F,didt,neib,dneib,possTest,dist,fixneib,fixdneib,possneib,possdneib)
    L = zeros(length(possTest),1); 
    disp(size(F))
    for node = 1:length(possTest)
        nnode = neib(node,:); nnode = nnode(nnode~=0);
        dnode = dneib(node,:); dnode = dnode(dnode~=0);
        fixn = fixneib(node,:); fixn = fixn(fixn~=0);
        fixd = fixdneib(node,:); fixd = fixd(fixd~=0);
        possn = possneib(node,:); possn = possn(possn~=0);
        possd = possdneib(node,:); possd = possd(possd~=0);
        indn = []; indd = [];
        if isempty(possn) == 1
            indn = [];
        else
            for i = 1:length(possn)
                indn(i) = find(possTest == possn(i));
            end
        end
        if isempty(possd) == 1
            indd = [];
        else
            for i = 1:length(possd)
                indd(i) = find(possTest == possd(i));
            end
        end
%         L(node) = ((0.5.*didt(node)'*I(node)) + (tau*0.5.*didt(nnode)'*I(nnode)) + ((tau^2)*0.5.*didt(dnode)'*I(dnode)))*...
%             (1 - M(1)*R(node) - G(1)*F(node) - M(2)*sum(R(nnode)) ...
%             - G(2)*sum(F(nnode)) - M(3)*sum(R(dnode)) - G(3)*sum(F(dnode)));
        L(node) = ((didt(nnode)'*I(nnode))/dist + (didt(dnode)'*I(dnode))/(2*dist))*...
            (1 - M(1)*0 - G(1)*F(node) - M(2)*length(fixn) - G(2)*sum(F(indn)) - M(3)*length(fixd) - G(3)*sum(F(indd)));
    end
    L = -F*L;
end