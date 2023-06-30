one = false;

if (one)
    img_path = "Scitech Final Craters/lunar_032gs_3hD_MGB.png";

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
    disp(ang);
else
    csv = "Scitech Final Craters/2022_11_01_50mTorr_h3_1s_032gs_noacrylic_X-Slice.csv";

    data = readmatrix(csv);

    x = data(:,1);
    y = data(:,3);

%     [~, ind] = max(y);
%     
    ind = find(y, 12, "last");
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
    disp(ang);
end