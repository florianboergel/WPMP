% Function: Rotate_Matrix
%
%  This function generates a coordinate transformation matrix, based on the 
%  axis and angle of rotation.
%---------------------------------
% Keywords:
%------------
% Input:
%-------------
%   angle               Rotation angle in degrees
%   axis                Axis of rotation ('x','y' or 'z')
%------------
% Output:
%-------------
%   A                   Transformation matrix
%------------
% Created: 
% Marijn van Dooren, June 2014
% (c) Universität Oldenburg
% ----------------------------------
function [A]= Rotate_Matrix(angle,axis)

theta=angle/180*pi;

switch axis
    case 'x'
        A=[1 0 0 ; 0 cos(theta) sin(theta) ; 0 -sin(theta) cos(theta)];
    case 'y'
        A=[cos(theta) 0 -sin(theta) ; 0 1 0 ; sin(theta) 0 cos(theta)];
    case 'z'   
        A=[cos(theta) sin(theta) 0 ; -sin(theta) cos(theta) 0 ; 0 0 1];
    otherwise
        disp('Specified axis does not exist!')
end