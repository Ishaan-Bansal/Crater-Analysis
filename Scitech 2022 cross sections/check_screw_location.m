% Load the image
image = imread('martian_86gs_10hD_MGS/Image1.jpg'); % Replace 'your_image.jpg' with the path to your image

% Display the image
imshow(image);

% Prompt the user to click on the image
fprintf('Click on the image to get pixel coordinates. Press Enter when done.\n');
[x, y] = ginput;

% Display the clicked pixel coordinates
fprintf('Clicked pixel coordinates:\n');
for i = 1:numel(x)
    fprintf('Point %d: x = %.2f, y = %.2f\n', i, x(i), y(i));
end
