function [b,h,d,Iz] = dim_sample(  group_number )
% function [b,h,d,Iz] = dim_sample(  group_number )
% input :
% group_number = A1,A3,B1,B2,B3,C1,C2,C3,D1,D2,D3,D4,E1,E
% outputs :
% b = width (mm)
% h = thickness (mm)
% d = distance between the support and the region of interest (ROI) (mm)
% Iz = planar second moment of area (or planar area moment of inertia)
% (mm^4)s

switch group_number
    case {'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10',...
            'C11','C12','C13','C14','C15','C16','C17','C18','C19','C20'}
        b = 50;
        h = 25;
        d = 25;
        
end
Iz = b*h^3/12;
end
