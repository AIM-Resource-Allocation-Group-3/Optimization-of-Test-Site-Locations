function [Adj,obj] = make_graph(m,n,xstart,ystart,xend,yend)
% Making graphs for overlayed on x-y axis
hx = 1/(m+1); hy = 1/(n+1);
x = zeros(ceil(1/hx*(xend-xstart))+1,1); y = zeros(floor(1/hy*(yend-ystart))+1,1);
for i = 1:length(x)
    x(i) = i*hx;
end
for j = 1:length(y)
    y(j) = (j-1)*hy;
end

sz = length(x)*length(y);
% Generate adjaceny matrix
Adj = diag(1*ones(1,sz-1),-1) + ...
    diag(1*ones(1,sz-1),1) + diag(1*ones(1,(length(x)-1)*length(y)),-length(y)) + ...
    diag(1*ones(1,(length(x)-1)*length(y)),length(y));

% If node is onl the boundary
for i = 1:length(x)-1
    Adj(i*length(y), i*length(y)+1) = 0;
    Adj(i*length(y)+1, i*length(y)) = 0;
end
% Make graph
G = graph(Adj);

% Adding spatial component to graphs
% Using spatialgraph2D: https://www.mathworks.com/matlabcentral/fileexchange/73630-spatialgraph2d

spatialx = []; spatialy = [];
count = 0; 
while xstart+count*hx < xend + hx
    spatialx = [spatialx, xstart*ones(1,floor(1/hy*(yend-ystart)+1))+count*hx];
    count = count + 1;
end
for i = 1:ceil(1/hx*(xend-xstart))+1
    spatialy = [spatialy, yend:-hy:ystart];
end
obj = spatialgraph2D(G,spatialx,spatialy);
pgon=polyshape(obj);
end

