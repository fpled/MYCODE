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
    case {'B1','B2','B3','B4','B5','B6','B7','B8','B9',...
          'B10','B11','B12','B13','B14','B15','B16','B17','B18',...
          'B19','B20','B21','B22','B23','B24','B25','B26','B27'}
        b = 50;
        h = 15;
        d = 20;           
end
Iz = b*h^3/12;
end
