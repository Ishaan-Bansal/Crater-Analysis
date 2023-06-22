% ------- v2 --------
% lunar_032gs_3hD_MGB
% lunar_032gs_10hD_MDS
% lunar_032gs_10hD_MGB
% lunar_032gs_10hD_TRI
% martian_032gs_10hD_MGB
% martian_86gs_10hD_MDS
% martian_86gs_10hD_MGB

savePath = "v2\martian_86gs_10hD_MGB Slices";

% Sort the image files based on frame number for video creation
sortedImageFiles = dir(fullfile(savePath, 'plot_*.jpg'));
frameNumbers = regexp({sortedImageFiles.name}, 'plot_Image(\d+)', 'tokens');
frameNumbers = cellfun(@(x) str2double(x{1}), frameNumbers);
[~, sortedIndices] = sort(frameNumbers);
sortedImageFiles = sortedImageFiles(sortedIndices);

video = VideoWriter('v2\exports\martian_86gs_10hD_MGB.mp4', 'MPEG-4'); % Create the video object
open(video); % Open the file for writing

% Loop through each sorted image file and write to the video
for i = 1:numel(sortedImageFiles)
    imagePath = fullfile(savePath, sortedImageFiles(i).name);
    image = imread(imagePath);
    writeVideo(video, image); % Write the image to the video file
end

close(video); % Close the file
