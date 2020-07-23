%% SIS model using graphs
% Using distance to find adjacent and secondary nodes
m = 5; n = 5;
beta = 0.6; %infection coeff.
gamma = 0.2;  %recovery coeff.
tau = 0.8;   %movement b/t nodes coeff.
% Make graph
[A] = make_graph_unif(m);
H = graph(A);
% Assign attributes to nodes
dist = 3;possTest = [];
for i = 1:m*n
    if mod(i,dist) == 0
        H.Nodes.Suceptible(i) = floor(100.*(rand(1,1)+1))';
        H.Nodes.Infected(i) = floor(10.*(rand(1,1)+1))';
    elseif mod(i,5) == 0
        possTest = [possTest,i];
    else
        H.Nodes.Suceptible(i) = 0;
        H.Nodes.Infected(i) = 0;
    end
end
H.Nodes.Population = H.Nodes.Suceptible + H.Nodes.Infected;
popNodes = find(H.Nodes.Suceptible);
fixedTest = [4,8,19];
disp(H.Nodes)

% For plotting S,I over time
S_tot = sum(H.Nodes.Suceptible);
I_tot = sum(H.Nodes.Infected);

% Find nearest pop node within dist radius (neib)
% Find pop nodes within 2*dist radius (dneib)
% Find nearest fixed test site (fixneib)
% Find fixed test site 2*dist away (fixdneib)
% Find nearest poss test sites (possneib)
% Find poss test site 2*dist away (possdneib)
for node = 1:length(possTest)
    ntemp = nearest(H,possTest(node),dist); dtemp = nearest(H,possTest(node),2*dist);
    fixtemp = intersect(ntemp,fixedTest); fixdtemp = intersect(dtemp,fixedTest);
    posstemp = intersect(ntemp,possTest); possdtemp = intersect(dtemp,possTest);
    ntemp = ntemp(find(H.Nodes.Population(ntemp)));
    fixneib(node,1:length(fixtemp)) = fixtemp; possneib(node,1:length(posstemp)) = posstemp;
    dtemp = dtemp(find(H.Nodes.Population(dtemp)));
    fixdneibdiff = setdiff(fixdtemp,fixtemp); possdneibdiff = setdiff(possdtemp,posstemp);
    fixdneib(node,1:length(fixdneibdiff)) = fixdneibdiff; possdneib(node,1:length(possdneibdiff)) = possdneibdiff;
    neib(node,1:length(ntemp)) = ntemp;
    dtempDiff = setdiff(dtemp,ntemp);
    dneib(node,1:length(dtempDiff)) = dtempDiff;
    ntemp = []; dtemp = [];
end

pop_possTest = union(popNodes',possTest);
s = setdiff(1:m*n,pop_possTest);



dt = 1;
final_time = 5;
numsteps = final_time/dt;
x_tot = zeros(numsteps,1); I_indtot = zeros(numsteps,1);
nonlinx_tot = zeros(numsteps,1);
I_neigh = zeros(m*n,1); S_neigh = zeros(m*n,1);
didt = zeros(m*n,1); dsdt = zeros(m*n,1); Ival_ind = [];

for t = 1:numsteps
    % Calculate Infected,Suceptible neighbors
    for i = 1:m*n
        for j = 1:m*n
            if A(i,j) ~= 0
                I_neigh(i) = I_neigh(i) + H.Nodes.Infected(j);
                S_neigh(i) = S_neigh(i) + H.Nodes.Suceptible(j);
            end  
        end
    end
    disp(I_neigh)
    % total pop at a node
    P = H.Nodes.Suceptible + H.Nodes.Infected + tau.*(I_neigh + S_neigh);

    %Solving the ODE
    dsdt(popNodes) = -beta.*(H.Nodes.Infected(popNodes).*H.Nodes.Suceptible(popNodes))./P(popNodes) + gamma.*H.Nodes.Infected(popNodes) - tau.*(I_neigh(popNodes).*H.Nodes.Suceptible(popNodes))./P(popNodes);
    didt(popNodes) = beta.*(H.Nodes.Infected(popNodes).*H.Nodes.Suceptible(popNodes))./P(popNodes) - gamma.*H.Nodes.Infected(popNodes) + tau.*(I_neigh(popNodes).*H.Nodes.Suceptible(popNodes))./P(popNodes);
    H.Nodes.Infected = H.Nodes.Infected + didt.*dt;
    H.Nodes.Suceptible = H.Nodes.Suceptible + dsdt.*dt;

    I_tot = [I_tot, sum(H.Nodes.Infected)];
    S_tot = [S_tot, sum(H.Nodes.Suceptible)];
    
    % nonlinear optimization
    [nonlinx,val,exitFlag,Output] = nonlinear_opt_site_graphs(H.Nodes.Infected,didt,neib,dneib,possTest,dist,fixneib,fixdneib,possneib,possdneib);
    indnonlinx = find(nonlinx);
    if size(indnonlinx,2) == 0
        nonlinx_tot(t) = 0;
    else
        nonlinx_tot(t) = indnonlinx;
    end
    Ival = [];
    for node = 1:length(possTest)
        nnode = neib(node,:); nnode = nnode(nnode~=0);
        dnode = dneib(node,:); dnode = dnode(dnode~=0);
        Ival = [Ival,((didt(nnode)'*H.Nodes.Infected(nnode))/dist + (didt(dnode)'*H.Nodes.Infected(dnode))/(2*dist))];
    end
    Ival_ind(t) = find(Ival == max(Ival)); 
%     indI = find(H.Nodes.Infected == max(H.Nodes.Infected));
%     I_indtot(t) = indI;
    
end

figure
plot(I_tot,'r')
hold on
plot(S_tot,'b')
hold off

% Look at nodes assigned for objective function and Infected number
disp(possTest(nonlinx_tot))
disp(possTest(Ival_ind))

