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

% Load the image
image = imread(['Scitech Final Craters/lunar_032gs_3hD_MGB.png']); % Replace 'your_image.jpg' with the path to your image

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
