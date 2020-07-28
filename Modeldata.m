% Model with data
% There are 3 possible locations considered
%   1. popPoints - population nodes that are the centroids of the counties
%   2. possTest - possible pop-up test sites (right now this is places
%      of worship until we get more data)
%   3. fixedTest - fixed testing locations (right now this is hospitals)
S = shaperead('Data/tl_2019_33_bg/tl_2019_33_bg.shp');
County = shaperead('Data/NH_Population_Density-shp/USA_Population_Density.shp');
PlacesOfWorship = shaperead('Data/New_Hampshire_Places_of_Worship-shp/New_Hampshire_Places_of_Worship.shp');
Hospitals = shaperead('Data/New_Hampshire_Hospitals-shp/9c19a9cd-6911-4892-8c7b-ea3639c50f76202049-1-j3z5s6.y65y.shp');

% Find centroids of census tracts
numTracts = struct2cell(S);
Cx = zeros(size(numTracts,2),1); Cy = zeros(size(numTracts,2),1);
polgon = polyshape();
for i = 1:size(numTracts,2)
    polgon(i) = polyshape(S(i).X,S(i).Y);
    [Cx(i),Cy(i)] = centroid(polgon(i));
end
popPoints = []; % query points
for i = 1:length(Cx)
    popPoints = [popPoints;Cx(i), Cy(i)];
end

% Find points for possTest and fixedTest
Cxworship = zeros(size(struct2cell(PlacesOfWorship),2),1);
Cyworship = zeros(size(struct2cell(PlacesOfWorship),2),1);
Cxhosp = zeros(size(struct2cell(Hospitals),2),1);
Cyhosp = zeros(size(struct2cell(Hospitals),2),1);
for i = 1:size(struct2cell(PlacesOfWorship),2)
    Cxworship(i) = PlacesOfWorship(i).X1; 
    Cyworship(i) = PlacesOfWorship(i).Y1;
end
for i = 1:size(struct2cell(Hospitals),2)
    Cxhosp(i) = Hospitals(i).X1; 
    Cyhosp(i) = Hospitals(i).Y1;
end
possTestPoints = []; fixedTestPoints = [];
for i = 1:length(Cxworship)
    possTestPoints = [possTestPoints;Cxworship(i), Cyworship(i)];
end
for i = 1:length(Cxhosp)
    fixedTestPoints = [fixedTestPoints;Cxhosp(i), Cyhosp(i)];
end
%% County(6)

[INPOLY, ONPOLY] = isinterior(polyshape(County(6).X,County(6).Y),possTestPoints(:,1),possTestPoints(:,2));
ind = find(INPOLY);
newpossTestPoints = [];
newpossTestPoints(:,1) = possTestPoints(ind,1);
newpossTestPoints(:,2) = possTestPoints(ind,2);
possTestPoints = newpossTestPoints;

figure
mapshow(County)
hold on
plot(possTestPoints(:,1),possTestPoints(:,2),'.r')
plot(popPoints(:,1),popPoints(:,2),'.g')
hold off
%%
% subpossTestPoints = [];
% for i = 1:size(possTestPoints,1)
%     if possTestPoints(i,1) >= -71.6 && possTestPoints(i,1) <= -71.4 && possTestPoints(i,2) >= 42.7 && possTestPoints(i,2) <= 43.05
%         subpossTestPoints(i,1) = possTestPoints(i,1);
%         subpossTestPoints(i,2) = possTestPoints(i,2);
%     end
% end
% 
% newpossTestPoints = [];
% newpossTestPoints(:,1) = subpossTestPoints(subpossTestPoints(:,1)~=0,1);
% newpossTestPoints(:,2) = subpossTestPoints(subpossTestPoints(:,2)~=0,2);
% possTestPoints = newpossTestPoints;

% subpopPoints = [];
% for i = 1:size(popPoints,1)
%     if popPoints(i,1) >= -71.6 && popPoints(i,1) <= -71.4 && popPoints(i,2) >= 42.7 && popPoints(i,2) <= 43.05
%         subpopPoints(i,1) = popPoints(i,1);
%         subpopPoints(i,2) = popPoints(i,2);
%     end
% end
% 
% newpopPoints = [];
% newpopPoints(:,1) = subpopPoints(subpopPoints(:,1)~=0,1);
% newpopPoints(:,2) = subpopPoints(subpopPoints(:,2)~=0,2);
% popPoints = newpopPoints;



%% Finding primary and secondary neighbors
dist = 0.01; %radius
% Find the nearest pop to the possTest
popneib = []; popdneib = [];
[popneibIdx, popD] = rangesearch(popPoints,possTestPoints,dist); %dist away
[popdneibIdx, popdD] = rangesearch(popPoints,possTestPoints,2*dist); %2*dist away
% Find the nearest fixed test site to the possTest
fixneib = []; fixdneib = [];
[fixneibIdx, fixD] = rangesearch(fixedTestPoints,possTestPoints,dist); %dist away
[fixdneibIdx, fixdD] = rangesearch(fixedTestPoints,possTestPoints,2*dist); %2*dist
% Find the nearest pop-up test site to the possTest
possneib = []; possdneib = [];
[possneibIdx, possD] = rangesearch(possTestPoints,possTestPoints,dist);
[possdneibIdx, possdD] = rangesearch(possTestPoints,possTestPoints,2*dist);

for i = 1:size(popneibIdx,1)
    % Population
    if size(cell2mat(popneibIdx(i)),2) == 0 %no neighbors dist away
        popneib(i,:) = 0;
    end
    if size(cell2mat(popdneibIdx(i)),2) == 0 %no neighbors 2*dist away
        popdneib(i,:) = 0;
    end
    
    popneibtemp = cell2mat(popneibIdx(i)); popdneibtemp = cell2mat(popdneibIdx(i));
    popDiff = setdiff(popdneibtemp,popneibtemp);
    popneib(i,1:length(popneibtemp)) = popneibtemp;
    popdneib(i,1:length(popDiff)) = popDiff;
    
    % Fixed testing locations
    if size(cell2mat(fixneibIdx(i)),2) == 0 %no neighbors dist away
        fixneib(i,:) = 0;
    end
    if size(cell2mat(fixdneibIdx(i)),2) == 0 %no neighbors 2*dist away
        fixdneib(i,:) = 0;
    end
    
    fixneibtemp = cell2mat(fixneibIdx(i)); fixdneibtemp = cell2mat(fixdneibIdx(i));
    fixDiff = setdiff(fixdneibtemp,fixneibtemp);
    fixneib(i,1:length(fixneibtemp)) = fixneibtemp;
    fixdneib(i,1:length(fixDiff)) = fixDiff;
    
    % Pop-up testing locations
    if size(cell2mat(possneibIdx(i)),2) == 0 %no neighbors dist away
        possneib(i,:) = 0;
    end
    if size(cell2mat(possdneibIdx(i)),2) == 0 %no neighbors 2*dist away
        possdneib(i,:) = 0;
    end
    
    possneibtemp = cell2mat(possneibIdx(i)); possdneibtemp = cell2mat(possdneibIdx(i));
    possneibtemp = possneibtemp(possneibtemp~=i);
    possDiff = setdiff(possdneibtemp,possneibtemp);
    possneib(i,1:length(possneibtemp)) = possneibtemp;
    possdneib(i,1:length(possDiff)) = possDiff;
end

%% For calculating I_neigh and S_neigh
ineib = []; idneib = [];
[ineibIdx, iD] = rangesearch(popPoints,popPoints,dist); %dist away
[idneibIdx, idD] = rangesearch(popPoints,popPoints,2*dist); %2*dist away
for i = 1:size(ineibIdx,1)
    if size(cell2mat(ineibIdx(i)),2) == 0 %no neighbors dist away
        ineib(i,:) = 0;
    end
    if size(cell2mat(idneibIdx(i)),2) == 0 %no neighbors 2*dist away
        idneib(i,:) = 0;
    end
    
    ineibtemp = cell2mat(ineibIdx(i)); idneibtemp = cell2mat(idneibIdx(i));
    ineibtemp = ineibtemp(ineibtemp~=i);
    iDiff = setdiff(idneibtemp,ineibtemp);
    ineib(i,1:length(ineibtemp)) = ineibtemp;
    idneib(i,1:length(iDiff)) = iDiff;   
end

%% Disease modeling
%compartments
sz = size(popPoints,1);
TotPop = zeros(sz,1);
% for i = 1:sz
%     TotPop(i) = S(i).TOTPOP10;
% end
TotPop = floor(100.*(rand(1,sz)+1))';
I = floor(10.*(rand(1,sz)+1))';
Sus = TotPop - I;
Q = zeros(size(TotPop));
%parameters
beta = 0.6; %infection coeff.
gamma = 0.2;  %recovery coeff.
tau = 0.8;   %movement b/t nodes coeff.

dt = 1;
final_time = 1;
numsteps = final_time/dt;
x_tot = zeros(numsteps,1); I_indtot = zeros(numsteps,1);
nonlinx_tot = [];
I_neigh = zeros(sz,1); S_neigh = zeros(sz,1);
didt = zeros(sz,1); dsdt = zeros(sz,1); dqdt = zeros(sz,1); Ival_ind = []; I_tot = []; S_tot = []; q_tot = [];
pop_ntempsite = zeros(sz,1); pop_dtempsite = zeros(sz,1);
for t = 1:numsteps
    
    for i = 1:sz
        it = ineib(i,:); it = it(it~=0);
        for j = 1:length(it)
            I_neigh(i) = I_neigh(i) + I(it(j));
            S_neigh(i) = S_neigh(i) + Sus(it(j));
        end
    end
    % total pop at a node
    P = TotPop + tau.*(I_neigh + S_neigh);

    %Solving the ODE
    dsdt = -beta.*(I.*Sus)./P + gamma.*I - tau.*(I_neigh.*Sus)./P + (1/7).*Q_tot;
    didt = beta.*(I.*Sus)./P - gamma.*I + tau.*(I_neigh.*Sus)./P - ...
        (lambda1.*TotPop.*logical(fixneibIdx)./fixD + lambda2.*TotPop.*logical(fixedneibId)./fixdD + ...
        lambda1.*TotPop.*logical(sum(pop_ntempsite,2))./fixD + lambda2.*TotPop.*logical(sum(pop_dtempsite,2)));
    dqdt = lambda1.*TotPop.*logical(fixneibIdx)./fixD + lambda2.*TotPop.*logical(fixedneibId)./fixdD + ...
        lambda1.*TotPop.*logical(sum(pop_ntempsite,2))./fixD + lambda2.*TotPop.*logical(sum(pop_dtempsite,2)) - (1/7).*Q_tot;
    I = I + didt.*dt;
    Sus = Sus + dsdt.*dt;
    Q = Q + dqdt.*dt;

    I_tot = [I_tot, sum(I)];
    S_tot = [S_tot, sum(Sus)];
    Q_tot = [Q_tot, sum(Q)];
    
    % nonlinear optimization
    [nonlinx,val,exitFlag,Output,pop,s] = nonlinear_opt_site_graphs(I,didt,popneib,popdneib,1:size(possTestPoints,1),dist,fixneib,fixdneib,possneib,possdneib);
    indnonlinx = find(nonlinx); 
    if size(indnonlinx,2) == 0
        nonlinx_tot(t) = 0;
    else
        nonlinx_tot(1:length(indnonlinx),t) = indnonlinx;
    end
    Ival = [];
    for node = 1:size(possTestPoints,1)
        nnode = popneib(node,:); nnode = nnode(nnode~=0);
        dnode = popdneib(node,:); dnode = dnode(dnode~=0);
        Ival = [Ival,((didt(nnode)'*I(nnode))/dist + (didt(dnode)'*I(dnode))/(2*dist))];
    end 
    Ival_ind(t) = find(Ival == max(Ival));
%     indI = find(H.Nodes.Infected == max(H.Nodes.Infected));
%     I_indtot(t) = indI;
    
for i = 1:length(popPoints)
   ntempsite = rangesearch(indnonlinx,popPoints(i),dist); dtemp = nearest(indnonlinx,popPoints(i),2*dist); 
   dtempsite = setdiff(dtemp,ntempsite);
   pop_ntempsite(i,1:length(ntempsite)) = ntempsite;
   pop_dtempsite(i,1:length(dtempsite)) = dtempsite;   
end



end

figure
plot(I_tot,'r')
hold on
plot(S_tot,'b')
hold off

% Look at nodes assigned for objective function and Infected number
disp((nonlinx_tot))
disp((Ival_ind))


%%
figure
mapshow(County)
hold on
plot(possTestPoints(:,1),possTestPoints(:,2),'.r')
plot(popPoints(:,1),popPoints(:,2),'.g')
plot(fixedTestPoints(:,1),fixedTestPoints(:,2),'.k')
plot(possTestPoints(indnonlinx,1),possTestPoints(indnonlinx,2),'+r')
plot(possTestPoints(Ival_ind,1),possTestPoints(Ival_ind,2),'+b')
hold off
