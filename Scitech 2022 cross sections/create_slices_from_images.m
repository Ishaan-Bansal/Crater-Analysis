% lunar_032gs_3hD_MGS
% lunar_032gs_3hD_TRI
% lunar_032gs_10hD_MGS
% martian_032gs_10hD_MGS
% martian_86gs_10hD_MDS
% martian_86gs_10hD_MGS

% ------- v2 --------
% lunar_032gs_3hD_MGB
% lunar_032gs_10hD_MDS
% lunar_032gs_10hD_MGB
% lunar_032gs_10hD_TRI
% martian_032gs_10hD_MGB
% martian_032gs_10hD_TRI
% martian_86gs_10hD_MDS
% martian_86gs_10hD_MGB

% Specify the folder path containing the images
folderPath = 'v2\martian_86gs_10hD_MGB/edge case';
savePath = "v2/martian_86gs_10hD_MGB Slices";
origin = [374.75 28.25];
height = 85;

% Get a list of all image files in the folder
imageFiles = dir(fullfile(folderPath, '*.jpg')); % Update the file extension as per your image format

% Loop through each image file
for i = 1:numel(imageFiles)
    % Read the image
    imagePath = fullfile(folderPath, imageFiles(i).name);
    image = imread(imagePath); 

    [i2 mask] = removeRedDot(image);
    
    % Call the function to extract the edge coordinates
    [xCoordinates, yCoordinates] = extractEdge(i2, height, mask);
    yCoordinates = - yCoordinates;
    
    xCoordinates = xCoordinates - origin(1);
    yCoordinates = yCoordinates - origin(2);

    data = [xCoordinates, yCoordinates];

    csvFileName = fullfile(savePath, [imageFiles(i).name(1:end-4), '.csv']);
    dlmwrite(csvFileName, data, ',');

    
    % Plot the x and y coordinates
    figure;
    plot(xCoordinates, yCoordinates, 'k.');
    axis equal;
    
    % Save the plot as an image
    outputImagePath = fullfile(savePath, ['plot_', imageFiles(i).name]);
    saveas(gcf, outputImagePath);
    
    % Close the figure
    close(gcf);
end
