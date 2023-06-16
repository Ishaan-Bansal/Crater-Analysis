% Read the CSV files
filename1 = "lunar_032gs_10hD_MGS Slices/Image213";
filename2 = "martian_032gs_10hD_MGS Slices/Image353";
filename3 = "lunar_032gs_3hD_TRI Slices/Image327";
filename4 = "lunar_032gs_3hD_MGS Slices/Image327";
data1 = csvread(filename1 + '.csv');
data2 = csvread(filename2 + '.csv');
data3 = csvread(filename3 + '.csv');
data4 = csvread(filename4 + '.csv');

xData1 = data1(:, 1);
yData1 = data1(:, 2);

xData2 = data2(:, 1);
yData2 = data2(:, 2);

xData3 = data3(:, 1);
yData3 = data3(:, 2);

xData4 = data4(:, 1);
yData4 = data4(:, 2);

% Plot the data
figure;
hold on;
plot(xData1, yData1, 'b.', 'DisplayName', filename1);
plot(xData2, yData2, 'r.', 'DisplayName', filename2);
plot(xData3, yData3, 'g.', 'DisplayName', filename3);
plot(xData4, yData4, 'm.', 'DisplayName', filename4);
xlabel('X');
ylabel('Y');
title('Overlay of X and Y Coordinates');
l = legend('Location', 'best');
set(l, 'Interpreter', 'none')
hold off;
