function [x,y] = extractEdge(I)
BW = im2bw(I);
dim = size(BW)
col = round(dim(2)/4);
row = min(find(BW(:,col)))
boundary = bwtraceboundary(BW,[row, col],'E');

x = boundary(:,2)
y = boundary(:,1)

mask = y < 200

x(mask) = []
y(mask) = []
end