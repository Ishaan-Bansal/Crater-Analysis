% lunar_032gs_3hD_MGS
% lunar_032gs_3hD_TRI
% lunar_032gs_10hD_MGS
% martian_032gs_10hD_MGS
% martian_86gs_10hD_MGS
% martian_86gs_10hD_MDS

% ------- v2 --------
% lunar_032gs_3hD_MGB
% lunar_032gs_10hD_MDS
% lunar_032gs_10hD_MGB
% lunar_032gs_10hD_TRI
% martian_032gs_10hD_MGB
% martian_032gs_10hD_TRI
% martian_86gs_10hD_MDS
% martian_86gs_10hD_MGB

savePath = "v2\martian_86gs_10hD_MDS Slices"

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
step = floor(numel(sortedImageFiles) / 4);
for i = 1:step:numel(sortedImageFiles)
    imagePath = fullfile(savePath, sortedImageFiles(i).name);
    data = csvread(imagePath);
    x = data(:, 1);
    y = data(:, 2);
    plot(x, y, '.', 'Color', colors(i, :));
end

imagePath = fullfile(savePath, sortedImageFiles(numel(sortedImageFiles)).name);
data = csvread(imagePath);
x = data(:, 1);
y = data(:, 2);
plot(x, y, '.', 'Color', colors(numel(sortedImageFiles), :));

x0 = 100;
y0 = -800;
width = 1600;
height = 200;
axis([x0 width y0 height]);

hold off;

% Create a legend
legendCell = cell(numel(sortedImageFiles), 1);
for i = 1:step:numel(sortedImageFiles)
    [~, imageName, ~] = fileparts(sortedImageFiles(i).name);
    if ~isempty(imageName) % Filter out empty filenames
        legendCell{i} = string(imageName); % Convert the filename to string scalar
    end
end
[~, imageName, ~] = fileparts(sortedImageFiles(numel(sortedImageFiles)).name);
if ~isempty(imageName) % Filter out empty filenames
    legendCell{numel(sortedImageFiles)} = string(imageName); % Convert the filename to string scalar
end
legendCell = legendCell(~cellfun('isempty', legendCell)); % Remove empty cells
legend(legendCell, 'Location', 'best');