function [P]=ellipsoid_par(init_center,R,t)
%%======================
% Get the final ellipsodal parameters using the solved affine transformation
% Input
%--------------
% inti_center     the initial spherical center
% transform.R     the solved affine matrix
% transform.t     the translation vector
%
% Output
%-------------
% P=(x_c, y_c, z_c, a, b, c, \alpha, \beta, \gamma)   the nine geometric ...
%                                                     parameters of the fitted ellipsoidal surface
%=======================


%center
center=t+R*init_center';

%three semi-axis length
B=R*R';
[V,D]=eig(B);
a=sqrt(D(1,1));
b=sqrt(D(2,2));
c=sqrt(D(3,3));

%rotation angles
rx1=atan2(-V(3,1),sqrt(V(1,1)^2+V(2,1)^2));
ry1=atan2(V(2,1),V(1,1));
rz1=atan2(V(3,2),V(3,3));

P=[center' a b c ry1 rx1 rz1];

end