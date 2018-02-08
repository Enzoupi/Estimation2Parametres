function [Y] = J(F)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
global Uobs
n=length(F)+1;
h=1/n;
% Le h provient de la norme L2 discrète
Y=h/2*sum( (direct(F)-Uobs).^2 );
end

