data = readtable("CSV Files\surface pressure.csv");

% titleText = sprintf('Chamber Pressure %.2f mTorr _ Nozzle Height %.2f hD', chamberPressureValues(i), nozzleHeightValues(j));
titleText = "Surface Pressure";

% Create separate plots for Depth, Diameter, and Volume against Surface
% pressure
figure;
hold on;
grid on;
plot(data.SurfacePressure_kPa_, data.Depth_mm_, 'ro', 'DisplayName', 'Depth');
errorbar(data.SurfacePressure_kPa_, data.Depth_mm_,data.DepthError_mm_, "ro");
hold off;

% Set plot title and labels for Depth
title("Depth");
xlabel('Surface Pressure (kPa)');
ylabel('Depth (mm)');

width = 1.1*max(data.SurfacePressure_kPa_);
height = 1.1*max(data.Depth_mm_);
axis([0 width 0 height]);

saveas(gcf, ['Surface Pressure (kPa) vs Depth.fig']);

% Create a new figure for Diameter
figure;
hold on;
grid on;
plot(data.SurfacePressure_kPa_, data.Diameter_mm_, 'mo', 'DisplayName', 'Diameter');
errorbar(data.SurfacePressure_kPa_, data.Diameter_mm_,data.DiameterError_mm_, "mo");
hold off;

% Set plot title and labels for Diameter
title('Diameter');
xlabel('Surface Pressure (kPa)');
ylabel('Diameter (mm)');

width = 1.1*max(data.SurfacePressure_kPa_);
height = 1.1*max(data.Diameter_mm_);
axis([0 width 0 height]);

saveas(gcf, ['Surface Pressure (kPa) vs Diameter.fig']);

% Create a new figure for Volume
figure;
hold on;
grid on;
plot(data.SurfacePressure_kPa_, data.Volume_mm_3_, 'bo', 'DisplayName', 'Volume');
errorbar(data.SurfacePressure_kPa_, data.Volume_mm_3_,data.VolumeError_mm_3_, "bo");
hold off;

% Set plot title and labels for Volume
title("Volume");
xlabel('Surface Pressure (kPa)');
ylabel('Volume (mm^3)');

width = 1.1*max(data.SurfacePressure_kPa_);
height = 1.1*max(data.Volume_mm_3_);
axis([0 width 0 height]);

saveas(gcf, ['Surface Pressure (kPa) vs Volume.fig']);

% Create a new figure for Ridge Height
figure;
hold on;
grid on;
plot(data.SurfacePressure_kPa_, data.RidgeHeight_mm_, 'ko', 'DisplayName', 'Ridge Height');
errorbar(data.SurfacePressure_kPa_, data.RidgeHeight_mm_,data.RidgeHeightError_mm_, "ko");
hold off;

% Set plot title and labels for Volume
title("Ridge Height");
xlabel('Surface Pressure (kPa)');
ylabel('Ridge Height (mm)');

width = 1.1*max(data.SurfacePressure_kPa_);
height_max = 1.1*max(data.RidgeHeight_mm_);
height_min = 0.9*min(data.RidgeHeight_mm_);
axis([0 width height_min height_max]);

saveas(gcf, ['Surface Pressure (kPa) vs Ridge Height.fig']);