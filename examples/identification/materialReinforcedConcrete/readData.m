function [displ,force] = readData(filename)

displ = [];
force = [];
fid = fopen(filename);
if fid>0
    for i=1:4
        fgetl(fid);
    end
    h = fscanf(fid,'%f',[2 Inf])';
    displ = h(:,1);
    force = h(:,2);

%     while ~feof(fid)
%         tline = fgetl(fid);
%         h = sscanf(tline,'%f');
%         if ~isempty(h)
%             displ = [displ; h(1)];
%             force = [force; h(2)];
%         end
%     end
    
    fclose(fid);
end
