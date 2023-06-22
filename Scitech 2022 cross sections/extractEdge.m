function [x,y] = extractEdge(I, thresh, extMask)
I = imadjustn(I);
BW = im2bw(I, 0.5);
BW = medfilt2(BW,[5 5]);

% imshow(BW);
dim = size(BW);
col = round(dim(2)/2);
row = min(find(BW(:,col)));
boundary = bwtraceboundary(BW,[row, col],'SE');

x = boundary(:,2);
y = boundary(:,1);

mask = y < thresh;

x(mask) = [];
y(mask) = [];


% Get the linear indices of coordinates falling within the mask
indices = sub2ind(size(extMask), y, x);

% Reassign the coordinates within the mask to NaN
x(ismember(indices, find(extMask))) = NaN;
y(ismember(indices, find(extMask))) = NaN;
% x = fillgaps(x, 80, 1);
% y = fillgaps(y);

x = fillmissing(x,"pchip");
y = fillmissing(y,"pchip");
end


% imhist()