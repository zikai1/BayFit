%%===================
% Draw ellipsoids by the 9 geometric parameters
% Input: elli=[XC YC ZC A B C PHI THETA PSI];
%        type 0: no outliers 1: outliers
%        num: the number of points
% Output: draw the ellipsoid

function [Pt]=drawEllipsoid(elli,type,num)
%DRAWELLIPSOID Draw a 3D ellipsoid.
%
%   drawEllipsoid(ELLI)
%   Displays a 3D ellipsoid on current axis. ELLI is given by:
%   [XC YC ZC A B C PHI THETA PSI],
%   where (XC, YC, ZC) is the ellipsoid center, A, B and C are the half
%   lengths of the ellipsoid main axes, and PHI THETA PSI are Euler angles
%   representing ellipsoid orientation, in degrees.
%
%   drawEllipsoid(..., 'drawEllipses', true)
%   Also displays the main 3D ellipses corresponding to XY, XZ and YZ
%   planes.
%
%
%   Example
%     figure; hold on;
%     drawEllipsoid([10 20 30   50 30 10   5 10 0]);
%     axis equal;
%
%     figure; hold on;
%     elli = [[10 20 30   50 30 10   5 10 0];
%     drawEllipsoid(elli, 'FaceColor', 'r', ...
%         'drawEllipses', true, 'EllipseColor', 'b', 'EllipseWidth', 3);
%     axis equal;
%
%   See also
%   spheres, drawSphere, inertiaEllipsoid, ellipsoid, drawTorus, drawCuboid 
%

% ------
% Author: David Legland
% e-mail: david.legland@grignon.inra.fr
% Created: 2011-03-12,    using Matlab 7.9.0.529 (R2009b)
% Copyright 2011 INRA - Cepia Software Platform.


%% Default values

% number of meridians
nPhi    = num;

% number of parallels
nTheta  = num;

%% Parse numerical inputs
xc  = elli(:,1);
yc  = elli(:,2);
zc  = elli(:,3);
a   = elli(:,4);
b   = elli(:,5);
c   = elli(:,6);
ellPhi   = elli(:,7);
ellTheta = elli(:,8);
ellPsi   = elli(:,9);


%% Coordinates computation

% convert unit basis to ellipsoid basis
sca     = createScaling3d(a, b, c);
rotZ    = createRotationOz(ellPhi);
rotY    = createRotationOy(ellTheta);
rotX    = createRotationOx(ellPsi);
tra     = createTranslation3d([xc yc zc]);

% concatenate transforms
trans   = tra * rotZ * rotY * rotX * sca;


%% parameterisation of ellipsoid

% spherical coordinates
theta   = linspace(0, pi, nTheta+1);
phi     = linspace(0, 2*pi, nPhi+1);

% convert to cartesian coordinates
sintheta = sin(theta);
x = cos(phi') * sintheta;
y = sin(phi') * sintheta;
z = ones(length(phi),1) * cos(theta);
NP=numel(x);
for i = 1:NP
    res = [x(i) y(i) z(i) 1] * trans';
    x(i) = res(1);
    y(i) = res(2);
    z(i) = res(3);
end


switch type
    case 0
        x=x(:);
        y=y(:);
        z=z(:);
        
        
        %%{
           %add outliers
           outIntensity = 100;
           x_out=randn(1,56)*outIntensity;
           y_out=randn(1,56)*outIntensity;
           z_out=randn(1,56)*outIntensity;
           
        %}
        
        %Show the input data
        figure
        hold on;
        plot3(x,y,z,'r.');
        plot3(x_out,y_out,z_out,'k+');
        axis equal;
        grid on;
        view(3);
        
        hold off;
        
        %Show the final fitted ellipsoidal surface
        figure
        hold on
        plot3(x,y,z,'r.');
        plot3(x_out,y_out,z_out,'k+');
        axis equal;
        grid on;
        view(3);
        x=[x;x_out'];
        y=[y;y_out'];
        z=[z;z_out'];
    case 1
        x=x(:);
        y=y(:);
        z=z(:);
    case 2
        surf(x,y,z);
        view(3);
end




Pt=[x y z];


end

