%% Looking at shape files
% S - counties of NH (polygons)
% A - POI of NH (represented as a point)
% B - Census Tract of NH (polygons)

S = shaperead('NH_Population_Density-shp/USA_Population_Density.shp');
% A = shaperead('new_hampshire_point_of_interest/new_hampshire_poi.shp');
B = shaperead('tl_2016_33_tract/tl_2016_33_tract.shp');

Hospitals = shaperead('New_Hampshire_Hospitals-shp/9c19a9cd-6911-4892-8c7b-ea3639c50f76202049-1-j3z5s6.y65y.shp');
PlacesOfWorship = shaperead('New_Hampshire_Places_of_Worship-shp/New_Hampshire_Places_of_Worship.shp');
Pop = shaperead('tl_2019_33_bg/tl_2019_33_bg.shp');

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

mapshow(Pop)
hold on
plot(Cxworship,Cyworship,'.b','MarkerSize',10)
plot(Cxhosp,Cyhosp,'.r','MarkerSize',10)

numTracts = struct2cell(Pop);
Cx = zeros(size(numTracts,2),1); Cy = zeros(size(numTracts,2),1);
polgon = polyshape();
for i = 1:size(numTracts,2)
    polgon(i) = polyshape(Pop(i).X,Pop(i).Y);
    [Cx(i),Cy(i)] = centroid(polgon(i));
end
plot(Cx,Cy,'.g','MarkerSize',10)
hold off
