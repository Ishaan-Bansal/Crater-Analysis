% lunar_032gs_3hD_MGS
% lunar_032gs_3hD_TRI
% lunar_032gs_10hD_MGS
% martian_032gs_10hD_MGS
% martian_86gs_10hD_MDS
% martian_86gs_10hD_MGS

% Specify the folder path containing the images
folderPath = 'v2\martian_032gs_10hD_MGB';
savePath = "v2\martian_032gs_10hD_MGB Slices"
origin = [302.75 37.25]
height = 82

% Get a list of all image files in the folder
imageFiles = dir(fullfile(folderPath, '*.jpg')); % Update the file extension as per your image format

% Loop through each image file
for i = 1:numel(imageFiles)
    % Read the image
    imagePath = fullfile(folderPath, imageFiles(i).name);
    image = imread(imagePath);
    disp(imagePath)
    
    % Call the function to extract the edge coordinates
    [xCoordinates, yCoordinates] = extractEdge(image, height);
    yCoordinates = - yCoordinates
    
    xCoordinates = xCoordinates - origin(1)
    yCoordinates = yCoordinates - origin(2)

    data = [xCoordinates, yCoordinates]

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
