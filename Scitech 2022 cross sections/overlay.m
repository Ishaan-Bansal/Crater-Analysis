% ------- v2 --------
% lunar_032gs_3hD_MGB
% lunar_032gs_10hD_MDS
% lunar_032gs_10hD_MGB
% lunar_032gs_10hD_TRI
% martian_032gs_10hD_MGB
% martian_032gs_10hD_TRI
% martian_86gs_10hD_MDS
% martian_86gs_10hD_MGB

% Read the CSV files
filename1 = "v2\lunar_032gs_3hD_MGB Slices/Image327";
filename2 = "v2\lunar_032gs_10hD_MDS Slices/Image213";
filename3 = "v2\lunar_032gs_10hD_MGB Slices/Image213";
filename4 = "v2\lunar_032gs_10hD_TRI Slices/Image213";
filename5 = "v2\martian_032gs_10hD_MGB Slices/Image263";
filename6 = "v2\martian_032gs_10hD_TRI Slices/Image263";
filename7 = "v2\martian_86gs_10hD_MGB Slices/Image217";
data1 = csvread(filename1 + '.csv');
data2 = csvread(filename2 + '.csv');
data3 = csvread(filename3 + '.csv');
data4 = csvread(filename4 + '.csv');
data5 = csvread(filename5 + '.csv');
data6 = csvread(filename6 + '.csv');
data7 = csvread(filename7 + '.csv');

xData1 = data1(:, 1);
yData1 = data1(:, 2);

xData2 = data2(:, 1);
yData2 = data2(:, 2);

xData3 = data3(:, 1);
yData3 = data3(:, 2);

xData4 = data4(:, 1);
yData4 = data4(:, 2);

xData5 = data5(:, 1);
yData5 = data5(:, 2);

xData6 = data6(:, 1);
yData6 = data6(:, 2);

xData7 = data7(:, 1);
yData7 = data7(:, 2);

% Plot the data
figure;
hold on;
plot(xData1, yData1, 'b.', 'DisplayName', filename1);
plot(xData2, yData2, 'r.', 'DisplayName', filename2);
plot(xData3, yData3, 'g.', 'DisplayName', filename3);
plot(xData4, yData4, 'm.', 'DisplayName', filename4);
plot(xData5, yData5, 'k.', 'DisplayName', filename5);
plot(xData6, yData6, 'y.', 'DisplayName', filename6);
plot(xData7, yData7, 'c.', 'DisplayName', filename7);
xlabel('X');
ylabel('Y');
title('Overlay of X and Y Coordinates');
l = legend('Location', 'best');
set(l, 'Interpreter', 'none')
axis equal
hold off;
