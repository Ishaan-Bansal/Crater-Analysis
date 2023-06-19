% lunar_032gs_10hD_MDS 
% martian_032gs_10hD_MGB

vid=VideoReader('scitech videos 19th june/martian_032gs_10hD_MGB.mp4');
 numFrames = vid.NumFrames;
 n=numFrames;
 for i = 1:2:n
 frames = read(vid,i);
 imwrite(frames,['v2\martian_032gs_10hD_MGB\','Image' int2str(i), '.jpg']);
%  im(i)=image(frames);
 end
