 vid=VideoReader('scitech videos/lunar_032gs_3hD_TRI.mp4');
 numFrames = vid.NumFrames;
 n=numFrames;
 for i = 1:2:n
 frames = read(vid,i);
 imwrite(frames,['lunar_032gs_3hD_TRI\','Image' int2str(i), '.jpg']);
%  im(i)=image(frames);
 end
