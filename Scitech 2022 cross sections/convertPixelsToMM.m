function [x, y] = convertPixelsToMM(xData1, yData1)
x = xData1 * 63 / 106;
y = yData1 * 122 / 199;
end