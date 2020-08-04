% Model with data
% There are 3 possible locations considered
%   1. popPoints - population nodes that are the centroids of the counties
%   2. possTestPoints - possible pop-up test sites (right now this is places
%      of worship until we get more data)
%   3. fixedTestPoints - fixed testing locations (right now this is hospitals)

S = shaperead('Data/tl_2016_33_tract/tl_2016_33_tract.shp'); %census tract shape file
County = shaperead('Data/NH_Population_Density-shp/USA_Population_Density.shp'); %county shape file
PlacesOfWorship = shaperead('Data/New_Hampshire_Places_of_Worship-shp/New_Hampshire_Places_of_Worship.shp'); %places of worship shape file
Hospitals = shaperead('Data/New_Hampshire_Hospitals-shp/9c19a9cd-6911-4892-8c7b-ea3639c50f76202049-1-j3z5s6.y65y.shp'); %hospitals shape file

% Population data for census tracts - not currently being used since there
% are entries with 0's
[NUM,TXT,RAW] = xlsread('Data/TractPopulation.csv');
geoidString = string(cell2mat(RAW(2:end,2)));
geoid = erase(geoidString,"14000US");
pop = cell2mat(RAW(2:end,3));

% To create new field for pop
C = num2cell(zeros(1,length(S)));
[S(:).POP] = deal(C{:});

% Assign pop to corresponding geoid
for i = 1:length(S)
    for j = 1:length(S)
        if strcmp(geoid(i),string(S(j).GEOID)) == 1
            S(j).POP = pop(i);
            break
        end
    end
end

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
%% If you only wanted to look at 1 county instead of the whole state...
% Specifically for County(6) - Hillsborough County
% If you want to look at the whole state, skip this cell

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

%% Finding primary and secondary neighbors for possible testing site
% primary neighbors - neighbors dist away from each possible testing site
% secondary neighbors - neighbors 2*dist away from each possible testing site
% NOTE: Using rangesearch makes everything cells - some of the code done is
%       because of this

dist = 0.01; %radius
% Find the nearest population node to the possTest
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
    % Finding population nodes dist and 2*dist away
    if size(cell2mat(popneibIdx(i)),2) == 0 %no neighbors dist away
        popneib(i,:) = 0;
    end
    if size(cell2mat(popdneibIdx(i)),2) == 0 %no neighbors 2*dist away
        popdneib(i,:) = 0;
    end
    
    popneibtemp = cell2mat(popneibIdx(i)); popdneibtemp = cell2mat(popdneibIdx(i)); %convert to array
    popDiff = setdiff(popdneibtemp,popneibtemp);%take set diff to only include neighbors in the region [dist,2*dist] away
    popneib(i,1:length(popneibtemp)) = popneibtemp;%primary pop neighbors
    popdneib(i,1:length(popDiff)) = popDiff;%secondary pop neighbors
    
    % Finding fixed testing locations dist and 2*dist away
    if size(cell2mat(fixneibIdx(i)),2) == 0 %no neighbors dist away
        fixneib(i,:) = 0;
    end
    if size(cell2mat(fixdneibIdx(i)),2) == 0 %no neighbors 2*dist away
        fixdneib(i,:) = 0;
    end
    
    fixneibtemp = cell2mat(fixneibIdx(i)); fixdneibtemp = cell2mat(fixdneibIdx(i));%convert to array
    fixDiff = setdiff(fixdneibtemp,fixneibtemp);%take set diff to only include neighbors in the region [dist,2*dist]
    fixneib(i,1:length(fixneibtemp)) = fixneibtemp;%primary fix test site neighbors
    fixdneib(i,1:length(fixDiff)) = fixDiff;%secondary fix test site neighbors
    
    % Finding other pop-up testing locations dist and 2*dist away
    if size(cell2mat(possneibIdx(i)),2) == 0 %no neighbors dist away
        possneib(i,:) = 0;
    end
    if size(cell2mat(possdneibIdx(i)),2) == 0 %no neighbors 2*dist away
        possdneib(i,:) = 0;
    end
    
    possneibtemp = cell2mat(possneibIdx(i)); possdneibtemp = cell2mat(possdneibIdx(i));%convert to array
    possneibtemp = possneibtemp(possneibtemp~=i);%don't include self poss test site
    possDiff = setdiff(possdneibtemp,possneibtemp);%take set diff to only include neighbors in the region [dist,2*dist]
    possneib(i,1:length(possneibtemp)) = possneibtemp;%primary poss test site neighbors
    possdneib(i,1:length(possDiff)) = possDiff;%secondary poss test site neighbors
end

%% For calculating pop nodes dist and 2*dist away
% Used in I_neigh and S_neigh
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
    
    ineibtemp = cell2mat(ineibIdx(i)); idneibtemp = cell2mat(idneibIdx(i));%convert to array
    ineibtemp = ineibtemp(ineibtemp~=i);%don't include self pop node
    iDiff = setdiff(idneibtemp,ineibtemp);%take set diff to only include neighbors in the region [dist,2*dist]
    ineib(i,1:length(ineibtemp)) = ineibtemp;%primary pop node neighbors
    idneib(i,1:length(iDiff)) = iDiff;%secondary pop node neighbors
    
end

% Finding fixed testing site locations closest to the popNodes
nfixsite = []; dfixsite = [];
[nfixsiteIdx,nfd] = rangesearch(fixedTestPoints,popPoints,dist); 
[dfixsiteIdx,dfd] = rangesearch(fixedTestPoints,popPoints,2*dist);
for i = 1:size(nfixsiteIdx,1)
    if size(cell2mat(nfixsiteIdx(i)),2) == 0 %no neighbors dist away
        nfixsite(i,:) = 0;
    end
    if size(cell2mat(dfixsiteIdx(i)),2) == 0 %no neighbors 2*dist away
        dfixsite(i,:) = 0;
    end
    
    nfixtemp = cell2mat(nfixsiteIdx(i)); dfixtemp = cell2mat(dfixsiteIdx(i));%convert cell to array
    dtempsite = setdiff(dfixtemp,nfixtemp);%take set diff to only include neighbors in the region [dist,2*dist]
    nfixsite(i,1:length(nfixtemp)) = nfixtemp;%primary fixed test sites
    dfixsite(i,1:length(dtempsite)) = dtempsite;%secondary fixed test sites

    
end
%% Disease modeling
%compartments
sz = size(popPoints,1);
TotPop = zeros(sz,1);
% From earlier - need more accurate data for population
% for i = 1:sz
%     TotPop(i) = S(i).POP;
% end
TotPop = floor(100.*(rand(1,sz)+1))';%Assign random pop to each node
I = floor(1.*(rand(1,sz)+1))';%Assign random infected pop to each node
Q = zeros(size(TotPop));%Start with 0 people in quarantine
Sus = TotPop - I - Q;
%parameters
beta = 0.6; %infection coeff.
gamma = 0.2;  %recovery coeff.
tau = 0.2;   %movement b/t nodes coeff.
lambda1 = 0.002; %percent of people who go to fixed site
lambda2 = 0.001; %percent of people who go to pop-up site

dt = 1;%time step - 1 day
final_time = 2;
numsteps = final_time/dt;
x_tot = zeros(numsteps,1); I_indtot = zeros(numsteps,1);
nonlinx_tot = [];
I_neigh = zeros(sz,1); S_neigh = zeros(sz,1);
didt = zeros(sz,1); dsdt = zeros(sz,1); dqdt = zeros(sz,1); I_tot = []; S_tot = []; Q_tot = [];
pop_ntempsite = zeros(size(popPoints,1),5); pop_dtempsite = zeros(size(popPoints,1),5);%5 available pop-up sites to be assigned
for t = 1:numsteps
    % Find Infected and Susceptible neighbors
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
    dsdt = -beta.*(I.*Sus)./P + gamma.*I - tau.*(I_neigh.*Sus)./P + (1/7).*Q;
    didt = beta.*(I.*Sus)./P - gamma.*I + tau.*(I_neigh.*Sus)./P - ...
        (lambda1.*I.*sum(logical(nfixsite),2)./dist + lambda1.*I.*sum(logical(dfixsite),2)./(2*dist) + ...
        lambda2.*I.*sum(logical(pop_ntempsite),2)./dist + lambda2.*I.*sum(logical(pop_dtempsite),2)./(2*dist));
    dqdt = lambda1.*I.*sum(logical(nfixsite),2)./dist + lambda1.*I.*sum(logical(dfixsite),2)./(2*dist) + ...
        lambda2.*I.*sum(logical(pop_ntempsite),2)./dist + lambda2.*I.*sum(logical(pop_dtempsite),2)./(2*dist) - (1/7).*Q;
    I = I + didt.*dt;
    Sus = Sus + dsdt.*dt;
    Q = Q + dqdt.*dt;

    % For plotting
    I_tot = [I_tot, sum(I)];
    S_tot = [S_tot, sum(Sus)];
    Q_tot = [Q_tot, sum(Q)];
    
    % nonlinear optimization
    [nonlinx,val,exitFlag,Output,pop,s] = nonlinear_opt_site_graphs(I,didt,popneib,popdneib,1:size(possTestPoints,1),dist,fixneib,fixdneib,possneib,possdneib);
    % Keeping track of which sites got assigned
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
    disp(t)%keeping track of time
    % Finding assigned pop-up test sites closest to popNodes
    [ntempsite,ntempDist] = rangesearch(possTestPoints(indnonlinx',:),popPoints,dist); 
    [dtemp,dtempDist] = rangesearch(possTestPoints(indnonlinx',:),popPoints,2*dist);    

    for i = 1:length(ntempsite)
        if size(cell2mat(ntempsite(i)),2) == 0
            popntempsite(i,:) = 0;
        end
        if size(cell2mat(dtemp(i)),2) == 0
            pop_dtempsite(i,:) = 0;
        end
        popntempsite = cell2mat(ntempsite(i)); popdtempsite = cell2mat(dtemp(i));%convert cell to array
        dtempsite = setdiff(popdtempsite,popntempsite);%take set diff to only include neighbors in the region [dist,2*dist]
        pop_ntempsite(i,1:length(popntempsite)) = popntempsite;%primary assigned pop-up site neighbors
        pop_dtempsite(i,1:length(dtempsite)) = dtempsite;%secondary assigned pop-up site neighbors
    end


end

figure
plot(I_tot,'r')
hold on
plot(S_tot,'b')
hold off


% Look at nodes assigned for objective function 
disp((nonlinx_tot))

%%
figure 
f=worldmap([42.6 44],[-72.6 -70.4]);
geoshow(S,'FaceColor', [1 1 1],'DefaultEdgeColor', 'b') 
plotm(possTestPoints(indi,2),possTestPoints(indi,1),'r.');

markerSize = 50;
scatterm(possTestPoints(indi,2), possTestPoints(indi,1), markerSize, Ival(indi)','Filled');
colorbar

%% Animation
picWidth = 0.03;%how big you want the car
a = VideoWriter('2dayMovieSIQ.avi');%name it whatever you want
a.FrameRate = 1;%movie plays 1 frame every 1 second
a.Quality = 100;% range from [0,100]
open(a)
for j = 1:numsteps
    figure('Renderer', 'painters', 'Position', [10 10 900 600])%adjusts figure size
    mapshow(County)
    hold on
    plot(possTestPoints(:,1),possTestPoints(:,2),'.r','MarkerSize',8)
    plot(popPoints(:,1),popPoints(:,2),'.g','MarkerSize',8)
    plot(fixedTestPoints(:,1),fixedTestPoints(:,2),'.k','MarkerSize',8)
    for x = 1:size(nonlinx_tot,1)
        
        if nonlinx_tot(x,j) == 0
            break
        else
            ih = image(imread('car.png'));
            ih.XData = [possTestPoints(nonlinx_tot(x,j),1)-picWidth,possTestPoints(nonlinx_tot(x,j),1)+picWidth];
            ih.YData = [possTestPoints(nonlinx_tot(x,j),2)+picWidth,possTestPoints(nonlinx_tot(x,j),2)-picWidth];
        end
    end
    drawnow
    hold off
    M(j) = getframe;
    writeVideo(a,M(j));
end
close(a)
%% List of Assigned Test Sites for each day in excel worksheet
% Fields being recorded:
%   1. DAY
%   2. NAME
%   3. ADDRESS
%   4. CITY
%   5. STATE
%   6. ZIP
%   7. County
headers = {'Day','Name', 'Address', 'City', 'State', 'Zip', 'County'};
numAssign = size(nonlinx_tot,1);
nonlinx_tot = nonlinx_tot(:);
val = strings(length(nonlinx_tot),7);
dayval = []; count = 0;

for i = 1:length(nonlinx_tot)
    if mod(i,numAssign) == 1
        dayval = [dayval; ones(numAssign,1) + count];
        count = count + 1;
    end
    val(i,1) = PlacesOfWorship(nonlinx_tot(i)).NAME;
    val(i,2) = PlacesOfWorship(nonlinx_tot(i)).ADDRESS;
    val(i,3) = PlacesOfWorship(nonlinx_tot(i)).CITY;
    val(i,4) = PlacesOfWorship(nonlinx_tot(i)).STATE;
    val(i,5) = PlacesOfWorship(nonlinx_tot(i)).ZIP;
    val(i,6) = PlacesOfWorship(nonlinx_tot(i)).COUNTY;
end
%name it whatever you want
xlswrite('AssignedTestSites.xls',headers,1,'A1')
xlswrite('AssignedTestSites.xls',dayval,1,'A2')
xlswrite('AssignedTestSites.xls',val,1,'B2') 
    
    
    
    
