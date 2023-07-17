% lunar_032gs_10hD_MDS 
% martian_032gs_10hD_MGB

% ------- v2 --------
% lunar_032gs_3hD_MGB
% lunar_032gs_10hD_MDS
% lunar_032gs_10hD_MGB
% lunar_032gs_10hD_TRI
% martian_032gs_10hD_MGB
% martian_032gs_10hD_TRI
% martian_86gs_10hD_MDS
% martian_86gs_10hD_MGB

vid=VideoReader('scitech videos 19th june/lunar_032gs_3hD_MGB.mp4');
 numFrames = vid.NumFrames;
 n=numFrames;
 for i = 1:2:n
 frames = read(vid,i);
 imwrite(frames,['v2\lunar_032gs_3hD_MGB\','Image' int2str(i), '.jpg']);
%  im(i)=image(frames);
 end
