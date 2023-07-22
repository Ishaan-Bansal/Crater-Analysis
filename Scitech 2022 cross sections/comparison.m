% ------- v2 --------
% lunar_032gs_3hD_MGB
% lunar_032gs_10hD_MDS
% lunar_032gs_10hD_MGB
% lunar_032gs_10hD_TRI
% martian_032gs_10hD_MGB
% martian_032gs_10hD_TRI
% martian_86gs_10hD_MDS
% martian_86gs_10hD_MGB

% https://uofi.box.com/s/pses8c44c54httvdzzirkp7zefxzzvc9

% Read the CSV files
filename1 = "Scitech Final Craters\adjust.png";
filename1ab = "v2\martian_86gs_10hD_MGB Slices/Image351";
filename1b = "v2\martian_86gs_10hD_MGB Slices/Image251";
filename1c = "v2\martian_86gs_10hD_MGB Slices/Image31";
filename2 = "lab tests/03_10_2023_6Torr_10hD_86gs_2_X-Slice";
% data1 = readmatrix(filename1 + '.csv');

origin = [374.00,45.00];
height = 140;
image = imread(filename1); 

[i2, mask] = removeRedDot(image);

% Call the function to extract the edge coordinates
[xCoordinates, yCoordinates] = extractEdge(i2, height, mask);
yCoordinates = - yCoordinates;

xCoordinates = xCoordinates - origin(1);
yCoordinates = yCoordinates - origin(2);
data1 = [xCoordinates, yCoordinates];
data1ab = readmatrix(filename1ab + '.csv');
data1b = readmatrix(filename1b + '.csv');
data1c = readmatrix(filename1c + '.csv');
data2 = readmatrix(filename2 + '.csv');

xData1 = data1(:, 1);
yData1 = data1(:, 2);

[xData1, yData1] = convertPixelsToMM_2(xData1, yData1);

yData1 = yData1 + 110;
xData1 = xData1 - 302;

xData1ab = data1ab(:, 1);
yData1ab = data1ab(:, 2);

[xData1ab, yData1ab] = convertPixelsToMM(xData1ab, yData1ab);
[xData1ab, index_sort] = sort(xData1ab);
yData1ab = yData1ab(index_sort);

[~, minIndex] = min(yData1ab);

yAdjust = -70;
xAdjust = xData1ab(minIndex);

yData1ab = yData1ab - yAdjust;
xData1ab = xData1ab - xAdjust;

xData1b = data1b(:, 1);
yData1b = data1b(:, 2);

[xData1b, yData1b] = convertPixelsToMM(xData1b, yData1b);

yData1b = yData1b - yAdjust;
xData1b = xData1b - xAdjust;

xData1c = data1c(:, 1);
yData1c = data1c(:, 2);

[xData1c, yData1c] = convertPixelsToMM(xData1c, yData1c);

yData1c = yData1c - yAdjust;
xData1c = xData1c - xAdjust;


xData2 = data2(:, 1);
yData2 = data2(:, 3);

[xData2, index_sort] = sort(xData2);
yData2 = yData2(index_sort);

[minValue, minIndex] = min(yData2);
yData2 = yData2 - yData2(1);
xData2 = xData2 - xData2(minIndex);

% Plot the data
figure;
hold on;
plot(xData1, yData1, 'k.', 'DisplayName', filename1);
pab = plot(xData1ab, yData1ab, '.', 'DisplayName', filename1ab);
pab.Color = [0.5 0.5 0.5];
pb = plot(xData1b, yData1b, '.', 'DisplayName', filename1b);
pb.Color = [0.75 0.75 0.75];
pc = plot(xData1c, yData1c, '.', 'DisplayName', filename1c);
pc.Color = [0.9 0.9 0.9];
plot(xData2, yData2, 'r.', 'DisplayName', filename2);
xlabel('X');
ylabel('Y');
clear title;
title('Comparison of UIUC and MSFC data');
l = legend('Location', 'best');
set(l, 'Interpreter', 'none')
axis equal
hold off;
