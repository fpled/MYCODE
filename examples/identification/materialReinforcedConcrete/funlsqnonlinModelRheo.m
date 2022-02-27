function f = funlsqnonlinModelRheo(x,force_exp,epsilon,angle,varargin)
% function f = funlsqnonlinModelRheo(x,force_exp,epsilon,angle,varargin)

strain = epsilon(:)*cos(angle)*[1 -1]; 
stress = zeros(size(strain));              
for i=1:2
    stress(:,i) = solveModelRheo(x,strain(:,i),varargin{:}); % [MPa]
end
sigma = (stress(:,1)-stress(:,2))/(2*cos(angle)); % [MPa]

S = x(1); % equivalent rod section [m^2]
force = sigma*S; % [MN]
f = force - force_exp;

end