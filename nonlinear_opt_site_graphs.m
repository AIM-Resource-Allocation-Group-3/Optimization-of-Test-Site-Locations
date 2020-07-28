% Very slow and puts testing sites on boundaries
%INPUT LIST
    %currentI: vector of infected at every node
    %didt: vector of didt at every node
    %neib: matrix, each row is the neighboring population nodes < dist
    %dneib: matrix, each row is the neighboring population nodes in 
        %(dist,2dist)
    %possTest: vector of possible nodes where a popup site can be placed
    %dist: scalar, measure to determine neighbors
    %fixneib:matrix, each row is the nodes where fixed test sites are 
         %<dist from possible test site location
    %fixdneib: matrix, each row is the nodes where fixed test sites are in
        %(dist, 2*dist) from possible test site location
    %possneib: matrix, each row is the nodes where possible test site
        %locations are <dist from possible test site location
    %possdneib: matrix, each row is the nodes where psosible test site 
        %locations are in (dist,2dist) from possible test site location
function[x,fval,exitFlag,Output,pop,scores] = nonlinear_opt_site_graphs(currentI,didt,neib,dneib,possTest,dist,fixneib,fixdneib,possneib,possdneib)
I = currentI; 
C = 5; %number of pop-up testing sites
numPossTest = length(possTest);

syms M(b)
M(b) = piecewise(b==1,3,b==2,2,b==3,1); %Assigns penalization value for close fixed sites
syms G(b)
G(b) = piecewise(b==1,3,b==2,2,b==3,1); %Assigns penalization value for clsoe mobile sites 

M = [M(1), M(2), M(3)];
G = [G(1), G(2), G(3)];

Aeq = []; %empty, no equality constraint
beq = []; %empty, no equality constraint

%Inequality constraints
A = ones(1,numPossTest);
b = C; %maximum number of test sites

%Upper and lower bounds for the solution
lb = zeros(numPossTest,1);
ub = ones(numPossTest,1);

nonlcon = []; %no nonlinear constraints
intcon = 1:numPossTest; %setting all solution variables as integers
initPop = zeros(100,numPossTest); %setting initial number of guesses to 100

fun = @(F) objfun(I,M,G,F,didt,neib,dneib,possTest,dist,fixneib,fixdneib,possneib,possdneib);
% options for ga() solver (go to Input->options for full list): https://www.mathworks.com/help/gads/ga.html
    % PlotFcn - displays plot as the optimization is working
    % MaxStallGenerations - one possible stopping criteria (looks at avg relative change)
    % MaxGerations - number of iterations before ga stops
options = optimoptions('ga','PopulationSize',300,'PlotFcn',{@gaplotbestf,@gaplotstopping},...
    'FunctionTolerance',10e-3,'MaxStallGenerations',20,'MaxGenerations',50,'InitialPopulationMatrix',initPop);
options = optimoptions(options,'UseVectorized',true);
[x,fval,exitFlag,Output,pop,scores] = ga(fun,numPossTest,A,b,Aeq,beq,lb,ub,nonlcon,intcon,options);
end

function L = objfun(I,M,G,F,didt,neib,dneib,possTest,dist,fixneib,fixdneib,possneib,possdneib)
    L = zeros(length(possTest),1); 
    disp(size(F))
    for node = 1:length(possTest)
        nnode = neib(node,:); nnode = nnode(nnode~=0);                      %list of neighboring population nodes
        dnode = dneib(node,:); dnode = dnode(dnode~=0);                     %list of secondar neighbor population nodes
        fixn = fixneib(node,:); fixn = fixn(fixn~=0);                       %list of neighboring fixed sites 
        fixd = fixdneib(node,:); fixd = fixd(fixd~=0);                      %list of secondary fixed sights
        possn = possneib(node,:); possn = possn(possn~=0);                  %list of neighboring possible test sites
        possd = possdneib(node,:); possd = possd(possd~=0);                 %list of secondary possible test sites
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
        %OBJECTIVE FUNCTION
        L(node) = ((didt(nnode)'*I(nnode))/dist + (didt(dnode)'*I(dnode))/(2*dist))*...
            (1 - M(1)*0 - G(1)*F(node) - M(2)*length(fixn) - G(2)*sum(F(indn)) - M(3)*length(fixd) - G(3)*sum(F(indd)));
    end
    L = -F*L;
end