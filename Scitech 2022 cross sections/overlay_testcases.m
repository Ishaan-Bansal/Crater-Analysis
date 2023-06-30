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
filename1 = "v2\lunar_032gs_10hD_MDS Slices/Image31";
filename2 = "lab tests/11_01_2022_50mTorr_10hD_032gs_X-Slice";
data1 = readmatrix(filename1 + '.csv');
data2 = readmatrix(filename2 + '.csv');

xData1 = data1(:, 1);
yData1 = data1(:, 2);

[xData1, yData1] = convertPixelsToMM(xData1, yData1);

[xData1, index_sort] = sort(xData1);
yData1 = yData1(index_sort);

[minValue, minIndex] = min(yData1);

yData1 = yData1 - yData1(1);
xData1 = xData1 - xData1(minIndex);

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
plot(xData1, yData1, 'b.', 'DisplayName', filename1);
plot(xData2, yData2, 'r.', 'DisplayName', filename2);
xlabel('X');
ylabel('Y');
title('Overlay of X and Y Coordinates');
l = legend('Location', 'best');
set(l, 'Interpreter', 'none')
axis equal
hold off;
