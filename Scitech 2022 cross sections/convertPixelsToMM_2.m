function [x, y] = convertPixelsToMM_2(xData1, yData1)
x = xData1 * 19 / 26;
y = yData1 * 58.85 / 71;
end