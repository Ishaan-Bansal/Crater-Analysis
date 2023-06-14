% Specify the folder path containing the images
folderPath = 'Crater 1';
savePath = "Crater 1 Slices"

% Get a list of all image files in the folder
imageFiles = dir(fullfile(folderPath, '*.jpg')); % Update the file extension as per your image format

% Loop through each image file
for i = 1:numel(imageFiles)
    % Read the image
    imagePath = fullfile(folderPath, imageFiles(i).name);
    image = imread(imagePath);
    
    % Call the function to extract the edge coordinates
    [xCoordinates, yCoordinates] = extractEdge(image);
    yCoordinates = - yCoordinates

    data = [xCoordinates, yCoordinates]

    csvFileName = fullfile(savePath, [imageFiles(i).name(1:end-4), '.csv']);
    dlmwrite(csvFileName, data, ',');

    
    % Plot the x and y coordinates
    figure;
    plot(xCoordinates, yCoordinates);
    axis equal;
    
    % Save the plot as an image
    outputImagePath = fullfile(savePath, ['plot_', imageFiles(i).name]);
    saveas(gcf, outputImagePath);
    
    % Close the figure
    close(gcf);
end
