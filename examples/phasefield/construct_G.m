function [G] = construct_G(ma,mb)
% function [G] = construct_G(ma,mb)

% Construction of the diagonal
G = diag(diag(ma*mb'));
% Line 1
G = G + [ma(1); zeros(5,1)]*[mb(1); 0; 0; mb(4); mb(5); 0]';
G = G + [ma(4); zeros(5,1)]*[0; 2*mb(4); 0; mb(1); 0; mb(5)]';
G = G + [ma(5); zeros(5,1)]*[0; 0; 2*mb(5); 0; mb(1); mb(4)]';
% Line 2
G = G + [0; ma(2); zeros(4,1)]*[0; mb(2); 0; mb(4); 0; mb(6)]';
G = G + [0; ma(4); zeros(4,1)]*[0; 0; 0; mb(2); mb(6); 0]';
G = G + [0; ma(6); zeros(4,1)]*[0; 0; 2*mb(6); 0; mb(1); mb(4)]';
% Line 3
G = G + [0; 0; ma(3); zeros(3,1)]*[0; 0; mb(3); 0; mb(5); mb(6)]';
G = G + [0; 0; ma(5); zeros(3,1)]*[0; 0; 0; mb(6); mb(3); 0]';
G = G + [0; 0; ma(6); zeros(3,1)]*[0; 0; 0; mb(5); 0; mb(3)]';
% Line 4
G = G + [zeros(3,1); ma(1); 0; 0]*[0; 0; 0; mb(2); mb(6); 0]';
G = G + [zeros(3,1); ma(4); 0; 0]*[zeros(5,1); mb(6)]';
G = G + [zeros(3,1); ma(5); 0; 0]*[zeros(4,1); mb(4); mb(2)]';
% Line 5
G = G + [zeros(4,1); ma(1); 0]*[zeros(4,1); mb(3); 0]';
G = G + [zeros(4,1); ma(4); 0]*[zeros(5,1); mb(3)]';
G = G + [zeros(4,1); ma(5); 0]*[zeros(5,1); mb(6)]';
% Line 6
G = G + [zeros(5,1); ma(2)]*[zeros(5,1); mb(3)]';
%
% Symmetrization
G = (G + transpose(G))/2;