data = readtable("full for cross sections.csv");
disp(data)

% % --------- Control Variables: Chamber Pressure, Nozzle Height ------------
% 
% chamberPressureValues = unique(data.ChamberPressure_mTorr_);
% nozzleHeightValues = unique(data.NozzleHeight_h_D_);
% 
% % Create a cell array to store the separate tables
% splitTables = cell(numel(chamberPressureValues), numel(nozzleHeightValues));
% 
% % Loop through the unique chamber pressure values
% for i = 1:numel(chamberPressureValues)
%     currentPressure = chamberPressureValues(i);
%     
%     % Loop through the unique nozzle height values
%     for j = 1:numel(nozzleHeightValues)
%         currentHeight = nozzleHeightValues(j);
%         
%         % Filter the table based on the current chamber pressure and nozzle height
%         splitTables{i, j} = data(data.ChamberPressure_mTorr_ == currentPressure & data.NozzleHeight_h_D_ == currentHeight, :);
%     end
% end
% 
% % Loop through the splitTables
% for i = 1:numel(chamberPressureValues)
%     for j = 1:numel(nozzleHeightValues)
%         % Get the current table
%         currentTable = splitTables{i, j};
% 
%         if size(currentTable, 1) < 3
%             continue
%         end
% 
%         folder_index = currentTable.Folder_Index;
%         labels = num2str(currentTable.ID);
%         disp("run")
%         disp(labels)
%         
%         plotCrossSections(folder_index, labels)
%         
%     end
% end


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

        folder_index = currentTable.Folder_Index;
        labels = num2str(currentTable.NozzleHeight_h_D_);
        disp("run")
        disp(labels)
        
        plotCrossSections(folder_index, labels)

    end
end
% 
% % --------- Control Variables: Nozzle Height, Flow Rate -------------------
% 
% nozzleHeightValues = unique(data.NozzleHeight_h_D_);
% flowRateValues = unique(data.FlowRate_gs_);
% 
% % Create a cell array to store the separate tables
% splitTables = cell(numel(nozzleHeightValues), numel(flowRateValues));
% 
% % Loop through the unique chamber pressure values
% for i = 1:numel(nozzleHeightValues)
%     currentHeight = nozzleHeightValues(i);
%     
%     % Loop through the unique nozzle height values
%     for j = 1:numel(flowRateValues)
%         currentRate = flowRateValues(j);
%         
%         % Filter the table based on the current chamber pressure and nozzle height
%         splitTables{i, j} = data(data.NozzleHeight_h_D_ == currentHeight & data.FlowRate_gs_ == currentRate, :);
%     end
% end
% 
% % Loop through the splitTables
% for i = 1:numel(nozzleHeightValues)
%     for j = 1:numel(flowRateValues)
%         % Get the current table
%         currentTable = splitTables{i, j};
% 
%         if size(currentTable, 1) < 3
%             continue
%         end
% 
%         
%         
%     end
% end
% 
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
%                 
%             end
%         end
%     end
% end
