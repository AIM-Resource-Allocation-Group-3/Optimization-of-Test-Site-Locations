%% SIS model using graphs
% Using distance to find adjacent and secondary nodes
m = 12; n = 12;
beta = 0.6; %infection coeff.
gamma = 0.2;  %recovery coeff.
tau = 0.8;   %movement b/t nodes coeff.
% Make graph
[A] = make_graph_unif(m);
H = graph(A);
% Assign attributes to nodes
dist = 3;
for i = 1:m*n
    if mod(i,dist) == 0
        H.Nodes.Suceptible(i) = floor(100.*(rand(1,1)+1))';
        H.Nodes.Infected(i) = floor(10.*(rand(1,1)+1))';
    else
        H.Nodes.Suceptible(i) = 0;
        H.Nodes.Infected(i) = 0;
    end
end
H.Nodes.Population = H.Nodes.Suceptible + H.Nodes.Infected;
disp(H.Nodes)

% For plotting S,I over time
S_tot = sum(H.Nodes.Suceptible);
I_tot = sum(H.Nodes.Infected);

% Find nearest pop node within dist radius (neib)
% Find pop nodes within 2*dist radius (dneib)
for node = 1:m*n
    ntemp = nearest(H,node,dist);
    dtemp = nearest(H,node,2*dist);
    ntemp = ntemp(find(H.Nodes.Population(ntemp)));
    dtemp = dtemp(find(H.Nodes.Population(dtemp)));
    neib(node,1:length(ntemp)) = ntemp;
    dtempDiff = setdiff(dtemp,ntemp);
    dneib(node,1:length(dtempDiff)) = dtempDiff;
    ntemp = []; dtemp = [];
end

dt = 1;
final_time = 10;
numsteps = final_time/dt;
x_tot = zeros(numsteps,1); I_indtot = zeros(numsteps,1);
nonlinx_tot = zeros(numsteps,1);
I_neigh = zeros(m*n,1); S_neigh = zeros(m*n,1);

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
    % total pop at a node
    P = H.Nodes.Suceptible + H.Nodes.Infected + tau.*(I_neigh + S_neigh);

    %Solving the ODE
    dsdt = -beta.*(H.Nodes.Infected.*H.Nodes.Suceptible)./P + gamma.*H.Nodes.Infected - tau.*(I_neigh.*H.Nodes.Suceptible)./P;
    didt = beta.*(H.Nodes.Infected.*H.Nodes.Suceptible)./P - gamma.*H.Nodes.Infected + tau.*(I_neigh.*H.Nodes.Suceptible)./P;
    H.Nodes.Infected = H.Nodes.Infected + didt.*dt;
    H.Nodes.Suceptible = H.Nodes.Suceptible + dsdt.*dt;
    disp(size(didt))

    I_tot = [I_tot, sum(H.Nodes.Infected)];
    S_tot = [S_tot, sum(H.Nodes.Suceptible)];
    
    % nonlinear optimization
    [nonlinx,val,exitFlag,Output] = nonlinear_opt_site_graphs(m,n,H.Nodes.Infected,neib,dneib);
    indnonlinx = find(nonlinx);
    if size(indnonlinx,2) == 0
        nonlinx_tot(t) = 0;
    else
        nonlinx_tot(t) = indnonlinx;
    end
    
    indI = find(H.Nodes.Infected == max(H.Nodes.Infected));
    I_indtot(t) = indI;
    
end

figure
plot(I_tot,'r')
hold on
plot(S_tot,'b')
hold off

