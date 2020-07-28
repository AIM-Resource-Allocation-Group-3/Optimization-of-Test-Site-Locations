%% SIS model using graphs
% Using distance to find adjacent and secondary nodes
m = 5; n = 5;       %Size of mesh
beta = 0.6;         %infection rate
gamma = 0.2;        %recovery rate
tau = 0.8;          %movement b/t nodes coeff.
lambda1 = 0;      %Testing rate for primary neighbor sites
lambda2 = 0;      %Testing rate for secondary neighbor sites


% Make graph
[A] = make_graph_unif(m);
H = graph(A);
% Assign attributes to nodes

dist = 3;possTest = [];                                         %List of nodes where popup sites can be
for i = 1:m*n
    H.Nodes.Quarentine(i) = 0;
    if mod(i,dist) == 0
        H.Nodes.Suceptible(i) = floor(100.*(rand(1,1)+1))';     %Number of suceptible people on each node
        H.Nodes.Infected(i) = floor(10.*(rand(1,1)+1))';        %Number of infected people on each node
    elseif mod(i,5) == 0
        possTest = [possTest,i];
    else
        H.Nodes.Suceptible(i) = 0;
        H.Nodes.Infected(i) = 0;
    end
end
H.Nodes.Population = H.Nodes.Suceptible + H.Nodes.Infected;
popNodes = find(H.Nodes.Suceptible);                            %Nodes with a population on them
fixedTest = [4,8,19];                                           %Nodes with a fixed test
disp(H.Nodes)

% For plotting S,I over time
S_tot = sum(H.Nodes.Suceptible);
I_tot = sum(H.Nodes.Infected);
Q_tot = 0;

%FOR POSSIBLE TEST SITES
% Find nearest pop node within dist radius (neib)
% Find pop nodes within 2*dist radius (dneib)
% Find nearest fixed test site (fixneib)
% Find fixed test site 2*dist away (fixdneib)
% Find nearest poss test sites (possneib)
% Find poss test site 2*dist away (possdneib)

%   neib: For a given node that is a possible test site finds other nodes with a
%           population within a distance dist
%   dneib: For a given node that is a possible test site finds other nodes with a
%           population further than dist, <2*dist
%   fixneib: For a given node that is a possible test site finds other nodes that are
%           a fixed testing site within a distance dist
%   fixdneib: For a given node that is a possible test site finds other nodes that
%           are a fixed testing site further than dist, <2*dist
%   possneib: For a given node that is a possible test site finds other nodes that
%           are  possible testing sites within a distance dist
%   possdneib: For a given node that is a possible test site finds other nodes that
%           are possible testing sites further than dist < 2*dist
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







%   pop_fixneib: For a given node that has a population finds other nodes that are a
%           fixed test site within a distance dist
%   pop_fixddneib: For a given node that has a population finds other nodes that are a
%           fixed test site further than dist, <2*dist
%   pop_possneib: For a given node that has a population finds other nodes that are
%           a possible testing site within a distance dist
%   pop_possdneib: For a given node that has a population finds other nodes that
%           are a possible testing site further than dist, <2*dist
for node = 1:length(popNodes)   
    ntemp = nearest(H,popNodes(node),dist); dtemp = nearest(H,popNodes(node),2*dist);
    fixtemp = intersect(ntemp,fixedTest); fixdtemp = intersect(dtemp,fixedTest);
    posstemp = intersect(ntemp,possTest); possdtemp = intersect(dtemp,possTest);
    ntemp = ntemp(find(H.Nodes.Population(ntemp))); %maps nearest nodes to nearest pop nodes
    pop_fixneib(node,1:length(fixtemp)) = fixtemp; pop_possneib(node,1:length(posstemp)) = posstemp;
    dtemp = dtemp(find(H.Nodes.Population(dtemp)));
    fixdneibdiff = setdiff(fixdtemp,fixtemp); possdneibdiff = setdiff(possdtemp,posstemp);
    pop_fixdneib(node,1:length(fixdneibdiff)) = fixdneibdiff; pop_possdneib(node,1:length(possdneibdiff)) = possdneibdiff;
    ntemp = []; dtemp = [];
end


pop_possTest = union(popNodes',possTest);   
s = setdiff(1:m*n,pop_possTest);




dt = 1;                                                     
final_time = 25;    
numsteps = final_time/dt;
x_tot = zeros(numsteps,1); I_indtot = zeros(numsteps,1);
nonlinx_tot = zeros(numsteps,1);
I_neigh = zeros(m*n,1); S_neigh = zeros(m*n,1);
didt = zeros(m*n,1); dsdt = zeros(m*n,1); Ival_ind = [];

%New Model 
dqdt = zeros(m*n,1)
pop_ntempsite = zeros(size(pop_possneib));
pop_dtempsite = zeros(size(pop_possneib));
for t = 1:numsteps
    % Calculate Infected,Suceptible neighbors
    for i = 1:m*n
        for j = 1:m*n
            if A(i,j) ~= 0
                I_neigh(i) = I_neigh(i) + H.Nodes.Infected(j);              %I_neigh(i) is the sum of all infected neighbors of node i
                S_neigh(i) = S_neigh(i) + H.Nodes.Suceptible(j);            %S_neigh(i) is the sum of all suceptible neighbors of node i
            end  
        end
    end
    disp(I_neigh)
    % total pop at a node
    P = H.Nodes.Suceptible + H.Nodes.Infected + tau.*(I_neigh + S_neigh);   %P is the total number of people on a node during 1 timestep

    %Solving the ODE
    dsdt(popNodes) = -beta.*(H.Nodes.Infected(popNodes).*H.Nodes.Suceptible(popNodes))./P(popNodes) ...                 %Infection rate
        + gamma.*H.Nodes.Infected(popNodes) - ...                                                                       %Recovery Rate
        tau.*(I_neigh(popNodes).*H.Nodes.Suceptible(popNodes))./P(popNodes)...                                          %Visitors being Infected on node
        + (1/7)*H.Nodes.Quarentine(popNodes);                                                                           %Return from quarentine
    
    didt(popNodes) =  beta.*(H.Nodes.Infected(popNodes).*H.Nodes.Suceptible(popNodes))./P(popNodes)...                                 
        - gamma.*H.Nodes.Infected(popNodes)...                                                                          
        + tau.*(I_neigh(popNodes).*H.Nodes.Suceptible(popNodes))./P(popNodes) ...                                       
        -(lambda1.*H.Nodes.Infected(popNodes).*sum(logical(pop_fixneib),2)./dist...                                     %Quarentine rate after visiting fixed site < dist away
        + lambda1.*H.Nodes.Infected(popNodes).*sum(logical(pop_fixdneib),2)./(2*dist) + ...                             %Quarentine rate after visiting fixed seit in (dist, 2dist) away
        lambda2.*H.Nodes.Infected(popNodes).*sum(pop_ntempsite,2)./dist ...                                             %Quarentine rate after visiting popup site < dist away
        + lambda2.*H.Nodes.Infected(popNodes).*sum(pop_dtempsite,2)./(2*dist));                                         %Quarentine rate after visitng popup site in (dist,2dist) away
    
    dqdt(popNodes) = - (1/7).*H.Nodes.Quarentine(popNodes)...                                                           
        +lambda1.*H.Nodes.Infected(popNodes).*sum(logical(pop_fixneib),2)./dist...                                      
        + lambda1.*H.Nodes.Infected(popNodes).*sum(logical(pop_fixdneib),2)./(2*dist) + ...
        lambda2.*H.Nodes.Infected(popNodes).*sum(pop_ntempsite,2)./dist ...
        + lambda2.*H.Nodes.Infected(popNodes).*sum(pop_dtempsite,2)./(2*dist)...
       
  
   
    %Explanation of dQdt
        % Takes populations and multiplies by number of fixed or popup
        % sites and scales by lambda1 or lambda2
        %no capacity, model assumes testing capacity <<< infected
        
     %Updating model   
    H.Nodes.Infected = H.Nodes.Infected + didt.*dt;
    H.Nodes.Suceptible = H.Nodes.Suceptible + dsdt.*dt;
    H.Nodes.Quarentine = H.Nodes.Quarentine + dqdt.*dt
    
    
    %Updating vector of totals
    I_tot = [I_tot, sum(H.Nodes.Infected)];
    S_tot = [S_tot, sum(H.Nodes.Suceptible)];
    Q_tot = [Q_tot, sum(H.Nodes.Quarentine)];
    
    % nonlinear optimization
    [nonlinx,val,exitFlag,Output] = nonlinear_opt_site_graphs(H.Nodes.Infected,didt,neib,...
                                                            dneib,possTest,dist,fixneib,...
                                                            fixdneib,possneib,possdneib);
    indnonlinx = find(nonlinx);

    
    %Can Delete?
    Ival = [];
    for node = 1:length(possTest)
        nnode = neib(node,:); nnode = nnode(nnode~=0);
        dnode = dneib(node,:); dnode = dnode(dnode~=0);
        Ival = [Ival,((didt(nnode)'*H.Nodes.Infected(nnode))/dist + (didt(dnode)'*H.Nodes.Infected(dnode))/(2*dist))];
    end


%pop_ntempsite: popup sites < dist from pop centers
%pop_dtempsite: popup sites < 2dist from popcenters
for i = 1:length(popNodes)
   ntemp = nearest(H,popNodes(node),dist); dtemp = nearest(H,popNodes(node),2*dist); 
   ntempsite = intersect(ntemp,indnonlinx); dtempsite = intersect(dtemp,indnonlinx);
   pop_ntempsite(node,1:length(ntempsite)) = ntempsite;
   pop_dtempsite(node,1:length(dtempsite)) = dtempsite;
   
   
end

    
end

figure
plot(I_tot,'r')
hold on
plot(S_tot,'b')
hold on
plot(Q_tot,'*')
xlabel('Days');
ylabel('Individuals');
title('\beta = 0.6, \gamma = 0.2  N = 5')
legend('Infected', 'Suceptible', 'Quarentine')

% Look at nodes assigned for objective function and Infected number
%disp(possTest(nonlinx_tot))
%disp(possTest(Ival_ind))
