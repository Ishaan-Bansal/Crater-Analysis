 vid=VideoReader('crater1.mp4');
 numFrames = vid.NumFrames;
 n=numFrames;
 for i = 1:2:n
 frames = read(vid,i);
 imwrite(frames,['Image' int2str(i), '.jpg']);
 im(i)=image(frames);
 end
