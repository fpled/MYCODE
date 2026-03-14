

[X,Y] = meshgrid(-5:0.5:5);
Z = Y.*sin(X) - X.*cos(Y);

FrameCount = 100; % number of frames
frame(1,FrameCount) = struct('cdata',[],'colormap',[]);
for i=1:FrameCount
    Zi = rand*Z;
    surf(X,Y,Zi)
    set(gca,'ZLim',[min(Z(:)), max(Z(:))])
    frame(i) = getframe(gcf);
end


filename = 'solution';
pathname = '.';
formats = {'avi','mp4'};
% FrameRate = 30; % rate of video playback in frames per second.
Duration = 10; % duration of the output file in seconds
FrameRate = FrameCount/Duration;
Quality = 100;

% Create movie file
mov = cell(1,length(formats));
for i=1:length(formats)
    if strcmp(formats{i},'avi')
        mov{i} = VideoWriter(fullfile(pathname,filename),'Motion JPEG AVI');
    elseif strcmp(formats{i},'mp4')
        mov{i} = VideoWriter(fullfile(pathname,filename),'MPEG-4');
    elseif strcmp(formats{i},'mj2')
        mov{i} = VideoWriter(fullfile(pathname,filename),'Motion JPEG 2000');
    end
    mov{i}.FrameRate = FrameRate;
    mov{i}.Quality = Quality;
    open(mov{i});
    writeVideo(mov{i},frame); % add the frames to the movie
    close(mov{i});
end
