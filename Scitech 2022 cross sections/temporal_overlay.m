savePath = "Crater 1 Slices"

% Sort the image files based on frame number for video creation
sortedImageFiles = dir(fullfile(savePath, '*.csv'));
frameNumbers = regexp({sortedImageFiles.name}, 'Image(\d+)', 'tokens');
frameNumbers = cellfun(@(x) str2double(x{1}), frameNumbers);
[~, sortedIndices] = sort(frameNumbers);
sortedImageFiles = sortedImageFiles(sortedIndices);

figure;
hold on;

colors = jet(numel(sortedImageFiles)); % Generate a colormap for the frames

% Loop through each sorted image file and plot the overlay
for i = 1:40:numel(sortedImageFiles)
    imagePath = fullfile(savePath, sortedImageFiles(i).name);
    data = csvread(imagePath);
    x = data(:, 1);
    y = data(:, 2);
    plot(x, y, '-', 'Color', colors(i, :));
end

x0 = 100;
y0 = -800;
width = 1600;
height = 200;
axis([x0 width y0 height]);

hold off;

% Create a legend
legendCell = cell(numel(sortedImageFiles), 1);
for i = 1:40:numel(sortedImageFiles)
    [~, imageName, ~] = fileparts(sortedImageFiles(i).name);
    if ~isempty(imageName) % Filter out empty filenames
        legendCell{i} = string(imageName); % Convert the filename to string scalar
    end
end
legendCell = legendCell(~cellfun('isempty', legendCell)); % Remove empty cells
legend(legendCell, 'Location', 'best');