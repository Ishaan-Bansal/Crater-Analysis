data = readtable("full.csv");

dir = "Figures/";

figure;
hold on;
grid on;

max_v = max(data.Volume_mm_3_);
max_error = max(data.VolumeError_mm_3_);

% red = [220,20,60];
% green = [0,255,0];
% 
% rve = data.VolumeError_mm_3_/max_error;
% ve = zeros(length(rve),3);
% for i = 1:length(rve)
%     if rve(i) < 0.5
%         ve(i,:) = red.*rve(i);
%     else 
%         ve(i,:) = green.*rve(i);
%     end
% end

colormapName = 'parula'; % Choose your desired colormap

rve = data.VolumeError_mm_3_ / max_error;

% Create a colormap indexed by your 'rve' values
cmap = colormap(gca, colormapName);
colormapIndex = round(rve * (size(cmap, 1) - 1)) + 1;
colormapIndex = colormapIndex / max(colormapIndex);



scatter3(data.FlowRate_gs_,data.ChamberPressure_mTorr_,data.NozzleHeight_h_D_,500*data.Volume_mm_3_/max_v,colormapIndex, "filled");
% scatter3(data.FlowRate_gs_,data.ChamberPressure_mTorr_,data.NozzleHeight_h_D_,100, ve);
view(40,35)

title('Crater Volume (mm^3)')
xlabel('Flow Rate (g/s)');
ylabel('Chamber Pressure (mTorr)');
zlabel('Nozzle Height (h/D)');
colorbar;