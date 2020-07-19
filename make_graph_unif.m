function [Adj] = make_graph_unif(m)
% Making graphs for uniform grid
% For a 3x3 graph, h = 1/(3+1) = 1/4
h = 1/(m+1);
x = zeros(1/h-1,1); y = zeros(1/h-1,1);
for i = 1:length(x)
    x(i) = i*h;
end
for j = 1:length(y)
    y(j) = (j-1)*h;
end

sz = length(x)*length(y);
% Generate adjaceny matrix
Adj = diag(1*ones(1,sz-1),-1) + ...
    diag(1*ones(1,sz-1),1) + diag(1*ones(1,length(x)*(length(y)-1)),-length(x)) + ...
    diag(1*ones(1,length(x)*(length(y)-1)),length(x));

% If node is onl the boundary
for i = 1:length(y)-1
    if mod(length(x),length(x)) == 0
        Adj(i*length(x), i*length(x)+1) = 0;
        Adj(i*length(x)+1, i*length(x)) = 0;
    end
end
% Make graph
G = graph(Adj);
plot(G)


% Adding spatial component to graphs
% Using spatialgraph2D: https://www.mathworks.com/matlabcentral/fileexchange/73630-spatialgraph2d

% spatialx = []; spatialy = [];
% for i = 1:1/h-1
%     spatialx = [spatialx, i*ones(1,1/h-1)];
%     spatialy = [spatialy, (1/h-1):-1:1];
% end
% obj = spatialgraph2D(G,spatialx,spatialy);
% pgon=polyshape(obj);
% 
% figure
% plot(obj)
% hold on
% plot(pgon)
% hold off
end

