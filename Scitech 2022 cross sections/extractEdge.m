function [x,y] = extractEdge(I, thresh)
I = imadjustn(I);
BW = im2bw(I, 0.5);
BW = medfilt2(BW,[5 5]);

imshow(BW)
dim = size(BW)
col = round(dim(2)/4);
row = max(find(BW(:,col)))
boundary = bwtraceboundary(BW,[row, col],'SE');

x = boundary(:,2)
y = boundary(:,1)

mask = y < thresh

x(mask) = []
y(mask) = []
end


% imhist()