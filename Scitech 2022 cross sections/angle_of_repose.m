one = true;

if (one)
    img_path = "Scitech Final Craters/martian_86gs_10hD_MGB.png";

    img = imread(img_path);
    
    [i2, mask] = removeRedDot(img);
    
    [x, y] = extractEdge(i2, 17, mask);
    y = - y;
    
    [val, ind] = min(y);
    
    mask = x > x(ind);
    
    x(mask) = [];
    y(mask) = [];
    
    
    coeffs = polyfit(x, y, 1);
    
    interpolatedX = x;
    interpolatedY = polyval(coeffs, interpolatedX);
    
    figure;
    plot(x, y, 'ro', 'MarkerSize', 10);
    figure;
    plot(interpolatedX, interpolatedY, 'b-', 'LineWidth', 3);
    ang = atan(coeffs(1));
    disp(rad2deg(ang));
else
    csv = "lab tests\03_10_2023_6Torr_10hD_86gs_2_X-Slice.csv";

    data = readmatrix(csv);

    x = data(:,1);
    y = data(:,3);
    
    ind = find(y, 21, "last"); % Change value according to ridge
    mask = x < x(ind(1));
    
    x(mask) = [];
    y(mask) = [];

    [~, ind] = min(y);
    
    mask = x > x(ind);
    
    x(mask) = [];
    y(mask) = [];
    
    coeffs = polyfit(x, y, 1);
    
    interpolatedX = x;
    interpolatedY = polyval(coeffs, interpolatedX);
    
    figure;
    plot(x, y, 'ro', 'MarkerSize', 10);
    figure;    
    plot(interpolatedX, interpolatedY, 'b-', 'LineWidth', 3);
    ang = atan(coeffs(1));
    disp(rad2deg(ang));
end