data = readtable("full.csv");

% --------- Control Variables: Chamber Pressure, Nozzle Height ------------

chamberPressureValues = unique(data.ChamberPressure_mTorr_);
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
        splitTables{i, j} = data(data.ChamberPressure_mTorr_ == currentPressure & data.NozzleHeight_h_D_ == currentHeight, :);
    end
end

% Loop through the splitTables
for i = 1:numel(chamberPressureValues)
    for j = 1:numel(nozzleHeightValues)
        % Get the current table
        currentTable = splitTables{i, j};

        if size(currentTable, 1) < 3
            continue
        end

        titleText = sprintf('Chamber Pressure %.2f mTorr _ Nozzle Height %.2f hD', chamberPressureValues(i), nozzleHeightValues(j));
        
        % Create separate plots for Depth, Diameter, and Volume against FlowRate
        figure;
        hold on;
        plot(currentTable.FlowRate_gs_, currentTable.Depth_mm_, 'ro', 'DisplayName', 'Depth');
        hold off;
        
        % Set plot title and labels for Depth
        title(titleText);
        xlabel('Flow Rate (gs)');
        ylabel('Depth (mm)');

        width = 1.1*max(currentTable.FlowRate_gs_);
        height = 1.1*max(currentTable.Depth_mm_);
        axis([0 width 0 height]);

        saveas(gcf, ['Flow Rate (gs) vs Depth; Control Variables - ' titleText '.svg']);
        
        % Create a new figure for Diameter
        figure;
        hold on;
        plot(currentTable.FlowRate_gs_, currentTable.Diameter_mm_, 'mo', 'DisplayName', 'Diameter');
        hold off;
        
        % Set plot title and labels for Diameter
        title(titleText);
        xlabel('Flow Rate (gs)');
        ylabel('Diameter (mm)');

        width = 1.1*max(currentTable.FlowRate_gs_);
        height = 1.1*max(currentTable.Diameter_mm_);
        axis([0 width 0 height]);

        saveas(gcf, ['Flow Rate (gs) vs Diameter; Control Variables - ' titleText '.svg']);
        
        % Create a new figure for Volume
        figure;
        hold on;
        plot(currentTable.FlowRate_gs_, currentTable.Volume_mm_3_, 'bo', 'DisplayName', 'Volume');
        hold off;
        
        % Set plot title and labels for Volume
        title(titleText);
        xlabel('Flow Rate (gs)');
        ylabel('Volume (mm^3)');

        width = 1.1*max(currentTable.FlowRate_gs_);
        height = 1.1*max(currentTable.Volume_mm_3_);
        axis([0 width 0 height]);

        saveas(gcf, ['Flow Rate (gs) vs Volume; Control Variables - ' titleText '.svg']);
        
    end
end

% --------- Control Variables: Chamber Pressure, Flow Rate ----------------

chamberPressureValues = unique(data.ChamberPressure_mTorr_);
flowRateValues = unique(data.FlowRate_gs_);

% Create a cell array to store the separate tables
splitTables = cell(numel(chamberPressureValues), numel(flowRateValues));

% Loop through the unique chamber pressure values
for i = 1:numel(chamberPressureValues)
    currentPressure = chamberPressureValues(i);
    
    % Loop through the unique nozzle height values
    for j = 1:numel(flowRateValues)
        currentRate = flowRateValues(j);
        
        % Filter the table based on the current chamber pressure and nozzle height
        splitTables{i, j} = data(data.ChamberPressure_mTorr_ == currentPressure & data.FlowRate_gs_ == currentRate, :);
    end
end

% Loop through the splitTables
for i = 1:numel(chamberPressureValues)
    for j = 1:numel(flowRateValues)
        % Get the current table
        currentTable = splitTables{i, j};

        if size(currentTable, 1) < 3
            continue
        end

        titleText = sprintf('Chamber Pressure %.2f mTorr, Flow Rate %.2f gs', chamberPressureValues(i), flowRateValues(j));
        
        % Create separate plots for Depth, Diameter, and Volume against FlowRate
        figure;
        hold on;
        plot(currentTable.NozzleHeight_h_D_, currentTable.Depth_mm_, 'ro', 'DisplayName', 'Depth');
        hold off;
        
        % Set plot title and labels for Depth
        title(titleText);
        xlabel('Nozzle Height (h/D)');
        ylabel('Depth (mm)');

        width = 1.1*max(currentTable.NozzleHeight_h_D_);
        height = 1.1*max(currentTable.Depth_mm_);
        axis([0 width 0 height]);

        saveas(gcf, ['Nozzle Height (hD) vs Depth; Control Variables - ' titleText '.svg']);
        
        % Create a new figure for Diameter
        figure;
        hold on;
        plot(currentTable.NozzleHeight_h_D_, currentTable.Diameter_mm_, 'mo', 'DisplayName', 'Diameter');
        hold off;
        
        % Set plot title and labels for Diameter
        title(titleText);
        xlabel('Nozzle Height (h/D)');
        ylabel('Diameter (mm)');

        width = 1.1*max(currentTable.NozzleHeight_h_D_);
        height = 1.1*max(currentTable.Diameter_mm_);
        axis([0 width 0 height]);

        saveas(gcf, ['Nozzle Height (hD) vs Diameter; Control Variables - ' titleText '.svg']);
        
        % Create a new figure for Volume
        figure;
        hold on;
        plot(currentTable.NozzleHeight_h_D_, currentTable.Volume_mm_3_, 'bo', 'DisplayName', 'Volume');
        hold off;
        
        % Set plot title and labels for Volume
        title(titleText);
        xlabel('Nozzle Height (h/D)');
        ylabel('Volume (mm^3)');

        width = 1.1*max(currentTable.NozzleHeight_h_D_);
        height = 1.1*max(currentTable.Volume_mm_3_);
        axis([0 width 0 height]);

        saveas(gcf, ['Nozzle Height (hD) vs Volume; Control Variables  - ' titleText '.svg']);

    end
end

% --------- Control Variables: Nozzle Height, Flow Rate -------------------

nozzleHeightValues = unique(data.NozzleHeight_h_D_);
flowRateValues = unique(data.FlowRate_gs_);

% Create a cell array to store the separate tables
splitTables = cell(numel(nozzleHeightValues), numel(flowRateValues));

% Loop through the unique chamber pressure values
for i = 1:numel(nozzleHeightValues)
    currentHeight = nozzleHeightValues(i);
    
    % Loop through the unique nozzle height values
    for j = 1:numel(flowRateValues)
        currentRate = flowRateValues(j);
        
        % Filter the table based on the current chamber pressure and nozzle height
        splitTables{i, j} = data(data.NozzleHeight_h_D_ == currentHeight & data.FlowRate_gs_ == currentRate, :);
    end
end

% Loop through the splitTables
for i = 1:numel(nozzleHeightValues)
    for j = 1:numel(flowRateValues)
        % Get the current table
        currentTable = splitTables{i, j};

        if size(currentTable, 1) < 3
            continue
        end

        titleText = sprintf('Nozzle Height %.2f hD, Flow Rate %.2f gs', nozzleHeightValues(i), flowRateValues(j));
        
        % Create separate plots for Depth, Diameter, and Volume against FlowRate
        figure;
        hold on;
        plot(currentTable.ChamberPressure_mTorr_, currentTable.Depth_mm_, 'ro', 'DisplayName', 'Depth');
        hold off;
        
        % Set plot title and labels for Depth
        title(titleText);
        xlabel('Chamber Pressure (mTorr)');
        ylabel('Depth (mm)');

        width = 1.1*max(currentTable.ChamberPressure_mTorr_);
        height = 1.1*max(currentTable.Depth_mm_);
        axis([0 width 0 height]);

        saveas(gcf, ['Chamber Pressure (mTorr) vs Depth; Control Variables - ' titleText '.svg']);
        
        % Create a new figure for Diameter
        figure;
        hold on;
        plot(currentTable.ChamberPressure_mTorr_, currentTable.Diameter_mm_, 'mo', 'DisplayName', 'Diameter');
        hold off;
        
        % Set plot title and labels for Diameter
        title(titleText);
        xlabel('Chamber Pressure (mTorr)');
        ylabel('Diameter (mm)');

        width = 1.1*max(currentTable.ChamberPressure_mTorr_);
        height = 1.1*max(currentTable.Diameter_mm_);
        axis([0 width 0 height]);

        saveas(gcf, ['Chamber Pressure (mTorr) vs Diameter; Control Variables - ' titleText '.svg']);
        
        % Create a new figure for Volume
        figure;
        hold on;
        plot(currentTable.ChamberPressure_mTorr_, currentTable.Volume_mm_3_, 'bo', 'DisplayName', 'Volume');
        hold off;
        
        % Set plot title and labels for Volume
        title(titleText);
        xlabel('Chamber Pressure (mTorr)');
        ylabel('Volume (mm^3)');

        width = 1.1*max(currentTable.ChamberPressure_mTorr_);
        height = 1.1*max(currentTable.Volume_mm_3_);
        axis([0 width 0 height]);

        saveas(gcf, ['Chamber Pressure (mTorr) vs Volume; Control Variables - ' titleText '.svg']);
        
    end
end

% % ---------------- Repeat Tests: Standard Deviation -----------------------
% 
% % Get unique values of Chamber Pressure, Nozzle Height, and Flow Rate
% chamberPressureValues = unique(data.ChamberPressure_mTorr_);
% nozzleHeightValues = unique(data.NozzleHeight_h_D_);
% flowRateValues = unique(data.FlowRate_gs_);
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
%             mask = (data.ChamberPressure_mTorr_ == currentChamberPressure) & ...
%                    (data.NozzleHeight_h_D_ == currentNozzleHeight) & ...
%                    (data.FlowRate_gs_ == currentFlowRate);
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
%                 titleText = sprintf('Chamber Pressure %.2f mTorr, Nozzle Height %.2f hD, Flow Rate %.2f gs', chamberPressureValues(i), nozzleHeightValues(j), flowRateValues(k));
%                 title(titleText);
%                 xlabel('Depth (mm)');
%                 ylabel('Probability Density');
%                 
%                 % Save the figure as .svg file with the specified naming convention
%                 saveas(gcf, ['Depth_NormalDistribution; Control Variables - ' titleText '.svg']);
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
%                 % Save the figure as .svg file with the specified naming convention
%                 saveas(gcf, ['Diameter_NormalDistribution; Control Variables - ' titleText '.svg']);
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
%                 % Save the figure as .svg file with the specified naming convention
%                 saveas(gcf, ['Volume_NormalDistribution; Control Variables - ' titleText '.svg']);
%             end
%         end
%     end
% end
% 
