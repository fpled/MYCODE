function f = funlsqnonlin(param,force_exp,epsilon,angle,varargin)
% function f = funlsqnonlin(param,force_exp,epsilon,angle,varargin)

strain = epsilon(:)*cos(angle)*[1 -1]; 
stress = zeros(size(strain));              
for i=1:2
    stress(:,i) = solveModelRheo(x,strain(:,i),varargin{:});
end
sigma = (stress(:,1)-stress(:,2))/(2*cos(angle));

S = x(1); % equivalent rod section [m^2]
force = sigma*S;
f = uy - uyexp;

end