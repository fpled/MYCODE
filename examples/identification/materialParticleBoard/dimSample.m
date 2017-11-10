function [b,h,d,Iz] = dimSample( sample )
% function [b,h,d,Iz] = dimSample( sample )
% input :
% sample = 'A1','A3','B1','B2','B3','C1','C2','C3','D1','D2','D3','D4','E1','E2'
% outputs :
% b = width (mm)
% h = thickness (mm)
% d = distance between the support and the region of interest (ROI) (mm)
% Iz = planar second moment of area (or planar area moment of inertia) (mm^4)

switch sample
    case {'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10',...
            'C11','C12','C13','C14','C15','C16','C17','C18','C19','C20'}
        b = 50;
        h = 25;
        d = 25;
        
end
Iz = b*h^3/12;
end
