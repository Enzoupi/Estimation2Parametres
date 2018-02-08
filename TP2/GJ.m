function [G] = GJ(F)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
U=direct(F);
G=adjoint(U);
end

