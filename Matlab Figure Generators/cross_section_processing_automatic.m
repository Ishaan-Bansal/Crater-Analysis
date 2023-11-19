data = readtable("CSV files/November 2023.csv");
disp(data)

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

        if size(currentTable, 1) < 3
            continue
        end

        folder_index = currentTable.Folder_Index;
        labels = cellstr(num2str(currentTable.FlowRate_g_s_));

        plotTitle = "Chamber Pressure " + num2str(currentTable{1,"ChamberPressure_Torr_"}) + " Torr, Nozzle Height " + num2str(currentTable{1,"NozzleHeight_h_D_"}) + " hD";
        plotCrossSections(folder_index, labels, "Flow Rate (g/s)", plotTitle)

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

    % Loop through the unique flow rate values
    for j = 1:numel(flowRateValues)
        currentRate = flowRateValues(j);

        % Filter the table based on the current chamber pressure and flow
        % rate
        splitTables{i, j} = data(data.ChamberPressure_Torr_ == currentPressure & data.FlowRate_g_s_ == currentRate, :);
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

        title = "Chamber Pressure " + num2str(currentTable{1,"ChamberPressure_Torr_"}) + " Torr, Flow Rate " + num2str(currentTable{1,"FlowRate_g_s_"}) + " gs";
        plotCrossSections(folder_index, labels, "Nozzle Height (h/D)", title)

    end
end

% --------- Control Variables: Nozzle Height, Flow Rate -------------------

nozzleHeightValues = unique(data.NozzleHeight_h_D_);
flowRateValues = unique(data.FlowRate_g_s_);

% Create a cell array to store the separate tables
splitTables = cell(numel(nozzleHeightValues), numel(flowRateValues));

% Loop through the unique nozzle height values
for i = 1:numel(nozzleHeightValues)
    currentHeight = nozzleHeightValues(i);

    % Loop through the unique flow rate values
    for j = 1:numel(flowRateValues)
        currentRate = flowRateValues(j);

        % Filter the table based on the current flow rate and nozzle height
        splitTables{i, j} = data(data.NozzleHeight_h_D_ == currentHeight & data.FlowRate_g_s_ == currentRate, :);
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

        folder_index = currentTable.Folder_Index;
        labels = num2str(currentTable.ChamberPressure_Torr_);

        title = "Nozzle Height " + num2str(currentTable{1,"NozzleHeight_h_D_"}) + " hD, Flow Rate " + num2str(currentTable{1,"FlowRate_g_s_"}) + " gs";
        plotCrossSections(folder_index, labels, "Chamber Pressure (Torr)", title)     

    end
end

% --------------------------- Repeat Tests --------------------------------

% Get unique values of Chamber Pressure, Nozzle Height, and Flow Rate
chamberPressureValues = unique(data.ChamberPressure_Torr_);
nozzleHeightValues = unique(data.NozzleHeight_h_D_);
flowRateValues = unique(data.FlowRate_g_s_);

% Loop through the unique combinations of Chamber Pressure, Nozzle Height, and Flow Rate
for i = 1:numel(chamberPressureValues)
    for j = 1:numel(nozzleHeightValues)
        for k = 1:numel(flowRateValues)
            % Get the current combination values
            currentChamberPressure = chamberPressureValues(i);
            currentNozzleHeight = nozzleHeightValues(j);
            currentFlowRate = flowRateValues(k);

            % Create a mask to filter the rows based on the current combination
            mask = (data.ChamberPressure_Torr_ == currentChamberPressure) & ...
                   (data.NozzleHeight_h_D_ == currentNozzleHeight) & ...
                   (data.FlowRate_g_s_ == currentFlowRate);

            % Get the subset of the table based on the mask
            currentTable = data(mask, :);

            % Check if the current table has less than 2 rows
            if size(currentTable, 1) < 2
                continue
            end
            folder_index = currentTable.Folder_Index;
            labels = currentTable.File;

            title = "Chamber Pressure " + num2str(currentTable{1,"ChamberPressure_Torr_"}) + " Torr, Nozzle Height " + num2str(currentTable{1,"NozzleHeight_h_D_"}) + " hD, Flow Rate " + num2str(currentTable{1,"FlowRate_g_s_"}) + " gs";
            plotCrossSections(folder_index, labels, "Filename", title)
        end
    end
end

% Individuals

for i = 1:numel(data.File)
    folder_index = data{i,"Folder_Index"};
    labels = data{i,"File"};

    title = "Chamber Pressure " + num2str(data{i,"ChamberPressure_Torr_"}) + " Torr, Nozzle Height " + num2str(data{i,"NozzleHeight_h_D_"}) + " hD, Flow Rate " + num2str(data{i,"FlowRate_g_s_"}) + " gs" + " ID " + data{i,"ID"};
    plotCrossSections(folder_index, labels, "Filename", title)
end

