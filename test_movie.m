

[X,Y] = meshgrid(-5:0.5:5);
Z = Y.*sin(X) - X.*cos(Y);
clear frame
for i=1:100
    Zi = rand*Z;
    surf(X,Y,Zi)
    set(gca,'zlim',[min(min(Z)),max(max(Z))])
    frame(i) = getframe(gcf);
end


filename = 'solution';
pathname = '.';
formats = {'avi','mp4'};
FrameRate = 30;
Quality = 100;

% Create movie file
mov = cell(1,length(formats));
for i=1:length(formats)
    if strcmp(formats{i},'avi')
        mov{i} = VideoWriter(fullfile(pathname,filename));
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
