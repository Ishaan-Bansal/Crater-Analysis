data = readtable("CSV Files/November 2023 III.csv");
% disp(data);

savefolder = "November 2023/";

% --------- Control Variables: Chamber Pressure, Nozzle Height ------------

chamberPressureValues = unique(data.ChamberPressure_Torr_);
nozzleHeightValues = unique(data.NozzleHeight_h_D_);

% Create a cell array to store the separate tables
splitTables = cell(numel(chamberPressureValues), numel(nozzleHeightValues));

% Loop through the unique chamber pressure values
for i = 1:numel(chamberPressureValues)
    currentPressure = chamberPressureValues(i);
    
    % Loop through the unique nozzle height values
    for j = 1:numel(nozzleHeightValues)
        currentHeight = nozzleHeightValues(j);
        
        % Filter the table based on the current chamber pressure and nozzle height
        splitTables{i, j} = data(data.ChamberPressure_Torr_ == currentPressure & data.NozzleHeight_h_D_ == currentHeight, :);
    end
end

% Loop through the splitTables
for i = 1:numel(chamberPressureValues)
    for j = 1:numel(nozzleHeightValues)
        % Get the current table
        currentTable = splitTables{i, j};

        if size(currentTable, 1) < 2
            continue
        end

        titleText = sprintf('Chamber Pressure %.2f Torr _ Nozzle Height %.2f hD', chamberPressureValues(i), nozzleHeightValues(j));
        
        % Create separate plots for Depth, Diameter, and Volume against FlowRate
        figure;
        hold on;
        grid on;
        plot(currentTable.FlowRate_g_s_, currentTable.Depth_mm_, 'ro', 'DisplayName', 'Depth');
        errorbar(currentTable.FlowRate_g_s_, currentTable.Depth_mm_,currentTable.DepthError_mm_, "ro");
        hold off;
        
        % Set plot title and labels for Depth
        title(titleText);
        xlabel('Flow Rate (g/s)');
        ylabel('Depth (mm)');

        width = 1.1*max(currentTable.FlowRate_g_s_);
        height = 1.1*max(currentTable.Depth_mm_);
        axis([0 width 0 height]);

        saveas(gcf, [savefolder 'Flow Rate (gs) vs Depth; Control Variables - ' titleText '.fig']);
        saveas(gcf, [savefolder 'Flow Rate (gs) vs Depth; Control Variables - ' titleText '.svg']);
        
        % Create a new figure for Diameter
        figure;
        hold on;
        grid on;
        plot(currentTable.FlowRate_g_s_, currentTable.Diameter_mm_, 'mo', 'DisplayName', 'Diameter');
        errorbar(currentTable.FlowRate_g_s_, currentTable.Diameter_mm_,currentTable.DiameterError_mm_, "mo");
        hold off;
        
        % Set plot title and labels for Diameter
        title(titleText);
        xlabel('Flow Rate (g/s)');
        ylabel('Diameter (mm)');

        width = 1.1*max(currentTable.FlowRate_g_s_);
        height = 1.1*max(currentTable.Diameter_mm_);
        axis([0 width 0 height]);

        saveas(gcf, [savefolder 'Flow Rate (g_s) vs Diameter; Control Variables - ' titleText '.fig']);
        saveas(gcf, [savefolder 'Flow Rate (g_s) vs Diameter; Control Variables - ' titleText '.svg']);
        
        % Create a new figure for Volume
        figure;
        hold on;
        grid on;
        plot(currentTable.FlowRate_g_s_, currentTable.Volume_mm_3_, 'bo', 'DisplayName', 'Volume');
        errorbar(currentTable.FlowRate_g_s_, currentTable.Volume_mm_3_,currentTable.VolumeError_mm_3_, "bo");
        hold off;
        
        % Set plot title and labels for Volume
        title(titleText);
        xlabel('Flow Rate (g/s)');
        ylabel('Volume (mm^3)');

        width = 1.1*max(currentTable.FlowRate_g_s_);
        height = 1.1*max(currentTable.Volume_mm_3_);
        axis([0 width 0 height]);

        saveas(gcf, [savefolder 'Flow Rate (g_s) vs Volume; Control Variables - ' titleText '.fig']);
        saveas(gcf, [savefolder 'Flow Rate (g_s) vs Volume; Control Variables - ' titleText '.svg']);

        % Create a new figure for Ridge Height
        figure;
        hold on;
        grid on;
        plot(currentTable.FlowRate_g_s_, currentTable.RidgeHeight_mm_, 'ko', 'DisplayName', 'Ridge Height');
        errorbar(currentTable.FlowRate_g_s_, currentTable.RidgeHeight_mm_,currentTable.RidgeHeightError_mm_, "ko");
        hold off;
        
        % Set plot title and labels for Volume
        title(titleText);
        xlabel('Flow Rate (g/s)');
        ylabel('Ridge Height (mm)');

        width = 1.1*max(currentTable.FlowRate_g_s_);
        height_max = 1.1*max(currentTable.RidgeHeight_mm_);
        height_min = 0.9*min(currentTable.RidgeHeight_mm_);
        axis([0 width height_min height_max]);

        saveas(gcf, [savefolder 'Flow Rate (g_s) vs Ridge Height; Control Variables - ' titleText '.fig']);
        saveas(gcf, [savefolder 'Flow Rate (g_s) vs Ridge Height; Control Variables - ' titleText '.svg']);
        
    end
end

% --------- Control Variables: Chamber Pressure, Flow Rate ----------------

chamberPressureValues = unique(data.ChamberPressure_Torr_);
flowRateValues = unique(data.FlowRate_g_s_);

% Create a cell array to store the separate tables
splitTables = cell(numel(chamberPressureValues), numel(flowRateValues));

% Loop through the unique chamber pressure values
for i = 1:numel(chamberPressureValues)
    currentPressure = chamberPressureValues(i);
    
    % Loop through the unique nozzle height values
    for j = 1:numel(flowRateValues)
        currentRate = flowRateValues(j);
        
        % Filter the table based on the current chamber pressure and nozzle height
        splitTables{i, j} = data(data.ChamberPressure_Torr_ == currentPressure & data.FlowRate_g_s_ == currentRate, :);
    end
end

% Loop through the splitTables
for i = 1:numel(chamberPressureValues)
    for j = 1:numel(flowRateValues)
        % Get the current table
        currentTable = splitTables{i, j};

        if size(currentTable, 1) < 2
            continue
        end

        titleText = sprintf('Chamber Pressure %.2f Torr, Flow Rate %.2f gs', chamberPressureValues(i), flowRateValues(j));
        
        % Create separate plots for Depth, Diameter, and Volume against FlowRate
        figure;
        hold on;
        grid on;
        plot(currentTable.NozzleHeight_h_D_, currentTable.Depth_mm_, 'ro', 'DisplayName', 'Depth');
        errorbar(currentTable.NozzleHeight_h_D_, currentTable.Depth_mm_,currentTable.DepthError_mm_, "ro");
        hold off;
        
        % Set plot title and labels for Depth
        title(titleText);
        xlabel('Nozzle Height (h/D)');
        ylabel('Depth (mm)');

        width = 1.1*max(currentTable.NozzleHeight_h_D_);
        height = 1.1*max(currentTable.Depth_mm_);
        axis([0 width 0 height]);

        saveas(gcf, [savefolder 'Nozzle Height (hD) vs Depth; Control Variables - ' titleText '.fig']);
        saveas(gcf, [savefolder 'Nozzle Height (hD) vs Depth; Control Variables - ' titleText '.svg']);
        
        % Create a new figure for Diameter
        figure;
        hold on;
        grid on;
        plot(currentTable.NozzleHeight_h_D_, currentTable.Diameter_mm_, 'mo', 'DisplayName', 'Diameter');
        errorbar(currentTable.NozzleHeight_h_D_, currentTable.Diameter_mm_,currentTable.DiameterError_mm_, "mo");
        hold off;
        
        % Set plot title and labels for Diameter
        title(titleText);
        xlabel('Nozzle Height (h/D)');
        ylabel('Diameter (mm)');

        width = 1.1*max(currentTable.NozzleHeight_h_D_);
        height = 1.1*max(currentTable.Diameter_mm_);
        axis([0 width 0 height]);

        saveas(gcf, [savefolder 'Nozzle Height (hD) vs Diameter; Control Variables - ' titleText '.fig']);
        saveas(gcf, [savefolder 'Nozzle Height (hD) vs Diameter; Control Variables - ' titleText '.svg']);
        
        % Create a new figure for Volume
        figure;
        hold on;
        grid on;
        plot(currentTable.NozzleHeight_h_D_, currentTable.Volume_mm_3_, 'bo', 'DisplayName', 'Volume');
        errorbar(currentTable.NozzleHeight_h_D_, currentTable.Volume_mm_3_,currentTable.VolumeError_mm_3_, "bo");
        hold off;
        
        % Set plot title and labels for Volume
        title(titleText);
        xlabel('Nozzle Height (h/D)');
        ylabel('Volume (mm^3)');

        width = 1.1*max(currentTable.NozzleHeight_h_D_);
        height = 1.1*max(currentTable.Volume_mm_3_);
        axis([0 width 0 height]);

        saveas(gcf, [savefolder 'Nozzle Height (hD) vs Volume; Control Variables  - ' titleText '.fig']);
        saveas(gcf, [savefolder 'Nozzle Height (hD) vs Volume; Control Variables  - ' titleText '.svg']);

         % Create a new figure for Ridge Height
        figure;
        hold on;
        grid on;
        plot(currentTable.NozzleHeight_h_D_, currentTable.RidgeHeight_mm_, 'ko', 'DisplayName', 'Ridge Height');
        errorbar(currentTable.NozzleHeight_h_D_, currentTable.RidgeHeight_mm_,currentTable.RidgeHeightError_mm_, "ko");
        hold off;
        
        % Set plot title and labels for Volume
        title(titleText);
        xlabel('Nozzle Height (h/D)');
        ylabel('Ridge Height (mm)');

        width = 1.1*max(currentTable.NozzleHeight_h_D_);
        height_max = 1.1*max(currentTable.RidgeHeight_mm_);
        height_min = 0.9*min(currentTable.RidgeHeight_mm_);
        axis([0 width height_min height_max]);

        saveas(gcf, [savefolder 'Nozzle Height (hD) vs Ridge Height; Control Variables - ' titleText '.fig']);
        saveas(gcf, [savefolder 'Nozzle Height (hD) vs Ridge Height; Control Variables - ' titleText '.svg']);
    end
end

% --------- Control Variables: Nozzle Height, Flow Rate -------------------

nozzleHeightValues = unique(data.NozzleHeight_h_D_);
flowRateValues = unique(data.FlowRate_g_s_);

% Create a cell array to store the separate tables
splitTables = cell(numel(nozzleHeightValues), numel(flowRateValues));

% Loop through the unique chamber pressure values
for i = 1:numel(nozzleHeightValues)
    currentHeight = nozzleHeightValues(i);
    
    % Loop through the unique nozzle height values
    for j = 1:numel(flowRateValues)
        currentRate = flowRateValues(j);
        
        % Filter the table based on the current chamber pressure and nozzle height
        splitTables{i, j} = data(data.NozzleHeight_h_D_ == currentHeight & data.FlowRate_g_s_ == currentRate, :);
    end
end

% Loop through the splitTables
for i = 1:numel(nozzleHeightValues)
    for j = 1:numel(flowRateValues)
        % Get the current table
        currentTable = splitTables{i, j};

        if size(currentTable, 1) < 2
            continue
        end

        titleText = sprintf('Nozzle Height %.2f hD, Flow Rate %.2f g_s', nozzleHeightValues(i), flowRateValues(j));
        
        % Create separate plots for Depth, Diameter, and Volume against FlowRate
        figure;
        hold on;
        grid on;
        plot(currentTable.ChamberPressure_Torr_, currentTable.Depth_mm_, 'ro', 'DisplayName', 'Depth');
        errorbar(currentTable.ChamberPressure_Torr_, currentTable.Depth_mm_,currentTable.DepthError_mm_, "ro");
        hold off;
        
        % Set plot title and labels for Depth
        title(titleText);
        xlabel('Chamber Pressure (Torr)');
        ylabel('Depth (mm)');

        width = 1.1*max(currentTable.ChamberPressure_Torr_);
        height = 1.1*max(currentTable.Depth_mm_);
        axis([0 width 0 height]);

        saveas(gcf, [savefolder 'Chamber Pressure (Torr) vs Depth; Control Variables - ' titleText '.fig']);
        saveas(gcf, [savefolder 'Chamber Pressure (Torr) vs Depth; Control Variables - ' titleText '.svg']);
        
        % Create a new figure for Diameter
        figure;
        hold on;
        grid on;
        plot(currentTable.ChamberPressure_Torr_, currentTable.Diameter_mm_, 'mo', 'DisplayName', 'Diameter');
        errorbar(currentTable.ChamberPressure_Torr_, currentTable.Diameter_mm_,currentTable.DiameterError_mm_, "mo");
        hold off;
        
        % Set plot title and labels for Diameter
        title(titleText);
        xlabel('Chamber Pressure (Torr)');
        ylabel('Diameter (mm)');

        width = 1.1*max(currentTable.ChamberPressure_Torr_);
        height = 1.1*max(currentTable.Diameter_mm_);
        axis([0 width 0 height]);

        saveas(gcf, [savefolder 'Chamber Pressure (Torr) vs Diameter; Control Variables - ' titleText '.fig']);
        saveas(gcf, [savefolder 'Chamber Pressure (Torr) vs Diameter; Control Variables - ' titleText '.svg']);
        
        % Create a new figure for Volume
        figure;
        hold on;
        grid on;
        plot(currentTable.ChamberPressure_Torr_, currentTable.Volume_mm_3_, 'bo', 'DisplayName', 'Volume');
        errorbar(currentTable.ChamberPressure_Torr_, currentTable.Volume_mm_3_,currentTable.VolumeError_mm_3_, "bo");
        hold off;
        
        % Set plot title and labels for Volume
        title(titleText);
        xlabel('Chamber Pressure (Torr)');
        ylabel('Volume (mm^3)');

        width = 1.1*max(currentTable.ChamberPressure_Torr_);
        height = 1.1*max(currentTable.Volume_mm_3_);
        axis([0 width 0 height]);

        saveas(gcf, [savefolder 'Chamber Pressure (Torr) vs Volume; Control Variables - ' titleText '.fig']);
        saveas(gcf, [savefolder 'Chamber Pressure (Torr) vs Volume; Control Variables - ' titleText '.svg']);

        % Create a new figure for Ridge Height
        figure;
        hold on;
        grid on;
        plot(currentTable.ChamberPressure_Torr_, currentTable.RidgeHeight_mm_, 'ko', 'DisplayName', 'Ridge Height');
        errorbar(currentTable.ChamberPressure_Torr_, currentTable.RidgeHeight_mm_,currentTable.RidgeHeightError_mm_, "ko");
        hold off;
        
        % Set plot title and labels for Volume
        title(titleText);
        xlabel('Chamber Pressure (Torr)');
        ylabel('Ridge Height (mm)');

        width = 1.1*max(currentTable.ChamberPressure_Torr_);
        height_max = 1.1*max(currentTable.RidgeHeight_mm_);
        height_min = 0.9*min(currentTable.RidgeHeight_mm_);
        axis([0 width height_min height_max]);

        saveas(gcf, [savefolder 'Chamber Pressure (Torr) vs Ridge Height; Control Variables - ' titleText '.fig']);        
        saveas(gcf, [savefolder 'Chamber Pressure (Torr) vs Ridge Height; Control Variables - ' titleText '.svg']);
    end
end

% % ---------------- Repeat Tests: Standard Deviation -----------------------
% 
% % Get unique values of Chamber Pressure, Nozzle Height, and Flow Rate
% chamberPressureValues = unique(data.ChamberPressure_Torr_);
% nozzleHeightValues = unique(data.NozzleHeight_h_D_);
% flowRateValues = unique(data.FlowRate_g_s_);
% 
% % Loop through the unique combinations of Chamber Pressure, Nozzle Height, and Flow Rate
% for i = 1:numel(chamberPressureValues)
%     for j = 1:numel(nozzleHeightValues)
%         for k = 1:numel(flowRateValues)
%             % Get the current combination values
%             currentChamberPressure = chamberPressureValues(i);
%             currentNozzleHeight = nozzleHeightValues(j);
%             currentFlowRate = flowRateValues(k);
%             
%             % Create a mask to filter the rows based on the current combination
%             mask = (data.ChamberPressure_Torr_ == currentChamberPressure) & ...
%                    (data.NozzleHeight_h_D_ == currentNozzleHeight) & ...
%                    (data.FlowRate_g_s_ == currentFlowRate);
%             
%             % Get the subset of the table based on the mask
%             currentTable = data(mask, :);
%             
%             % Check if the current table has more than 2 rows
%             if size(currentTable, 1) > 2
%                 % Create a new figure for Depth
%                 figure;
%                 pd = fitdist(currentTable.Depth_mm_, 'Normal');
%                 % Plot the PDF of Depth
%                 %histogram(currentTable.Depth_mm_, 'Normalization', 'pdf', 'FaceColor', 'b');
%                 scatter(currentTable.Depth_mm_, pdf(pd, currentTable.Depth_mm_), 'b');
%                 hold on;
%                 % Fit the normal distribution to Depth and plot the fitted curve
%                 x = linspace(min(currentTable.Depth_mm_), max(currentTable.Depth_mm_), 100);
%                 y = pdf(pd, x);
%                 plot(x, y, 'r', 'LineWidth', 2);
%                 hold off;
%                 
%                 % Set plot title and labels for Depth
%                 titleText = sprintf('Chamber Pressure %.2f Torr, Nozzle Height %.2f hD, Flow Rate %.2f g/s', chamberPressureValues(i), nozzleHeightValues(j), flowRateValues(k));
%                 title(titleText);
%                 xlabel('Depth (mm)');
%                 ylabel('Probability Density');
%                 
%                 % Save the figure as .fig file with the specified naming convention
%                 saveas(gcf, [savefolder 'Depth_NormalDistribution; Control Variables - ' titleText '.fig']);
%                 saveas(gcf, [savefolder 'Depth_NormalDistribution; Control Variables - ' titleText '.svg']);
%     
%                 % Create a new figure for Diameter
%                 figure;
%                 pd = fitdist(currentTable.Diameter_mm_, 'Normal');
%                 % Plot the PDF of Diameter
%                 %histogram(currentTable.Diameter_mm_, 'Normalization', 'pdf', 'FaceColor', 'g');
%                 scatter(currentTable.Diameter_mm_, pdf(pd, currentTable.Diameter_mm_), 'g');
%                 hold on;
%                 % Fit the normal distribution to Diameter and plot the fitted curve
%                 x = linspace(min(currentTable.Diameter_mm_), max(currentTable.Diameter_mm_), 100);
%                 y = pdf(pd, x);
%                 plot(x, y, 'r', 'LineWidth', 2);
%                 hold off;
%                 
%                 % Set plot title and labels for Diameter
%                 title(titleText);
%                 xlabel('Diameter (mm)');
%                 ylabel('Probability Density');
%                 
%                 % Save the figure as .fig file with the specified naming convention
%                 saveas(gcf, [savefolder 'Diameter_NormalDistribution; Control Variables - ' titleText '.fig']);
%                 saveas(gcf, [savefolder 'Diameter_NormalDistribution; Control Variables - ' titleText '.svg']);
%     
%                 % Create a new figure for Volume
%                 figure;
%                 pd = fitdist(currentTable.Volume_mm_3_, 'Normal');
%                 % Plot the PDF of Volume
%                 %histogram(currentTable.Volume_mm_3_, 'Normalization', 'pdf', 'FaceColor', 'm');
%                 scatter(currentTable.Volume_mm_3_, pdf(pd, currentTable.Volume_mm_3_), 'm');
%                 hold on;
%                 % Fit the normal distribution to Volume and plot the fitted curve
%                 x = linspace(min(currentTable.Volume_mm_3_), max(currentTable.Volume_mm_3_), 100);
%                 y = pdf(pd, x);
%                 plot(x, y, 'r', 'LineWidth', 2);
%                 hold off;
%                 
%                 % Set plot title and labels for Volume
%                 title(titleText);
%                 xlabel('Volume  (mm^3)');
%                 ylabel('Probability Density');
%                 
%                 % Save the figure as .fig file with the specified naming convention
%                 saveas(gcf, [savefolder 'Volume_NormalDistribution; Control Variables - ' titleText '.fig']);
%                 saveas(gcf, [savefolder 'Volume_NormalDistribution; Control Variables - ' titleText '.svg']);
%             end
%         end
%     end
% end
% 
