function [time,displ,force] = readDataCyclicCompression(filename)
% function [time,delta,force] = readDataCyclicCompression(filename)

data = load(filename);
time = data(:,1);
displ = -data(:,2);
force = -data(:,3);
displ = displ-displ(1);

% fid = fopen(filename,'r');
% if fid>0
%     data = fscanf(fid,'%f',[3 Inf])';
%     time = data(:,1);
%     displ = -data(:,2);
%     force = -data(:,3);
%     displ = displ-displ(1);
%     fclose(fid);
% end

end
