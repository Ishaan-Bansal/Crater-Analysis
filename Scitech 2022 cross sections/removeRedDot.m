function resultImage = removeRedDot(image)
    % Convert the image to the RGB color space
    rgbImage = im2double(image);
    
    % Extract the red, green, and blue channels
    redChannel = rgbImage(:, :, 1);
    greenChannel = rgbImage(:, :, 2);
    blueChannel = rgbImage(:, :, 3);
    
    % Threshold for red color detection (adjust as needed)
    redThreshold = 0.25;
    
    % Find the red pixels in the image
    redPixels = redChannel > redThreshold & ...
                greenChannel < redThreshold & ...
                blueChannel < redThreshold;
    
    % Create a mask of the red pixels
    mask = repmat(redPixels, [1, 1, 3]);
    
    % Replace the red pixels with the average of the surrounding pixels
    blurredImage = imgaussfilt(rgbImage, 30); % Adjust the standard deviation (5 in this case)
    resultImage = rgbImage .* (1 - mask) + blurredImage .* mask;
end
