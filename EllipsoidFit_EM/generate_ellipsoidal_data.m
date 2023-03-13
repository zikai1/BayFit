function [Samples1,L]=generate_ellipsoidal_data(ell_par,type,num)
%==========================================================================
% Generate ellipsoidal data points for fitting tests
% Method: First generate spherical points and then using the affine
% transformation to transform them as ellipsoidal points

% Input
%----------
% ell_par:  1x9  [Xc Yc Zc a b c alpha beta gamma]   ellipsoidal parameters;
% type:     0    without outliers; 
%           1    with outliers.
% num:      the number of data points.
%
% Output:
%----------
% Samples:  Nx3 array   the generated ellipsoidal points
% L:        3x3 affine matrix
%%=========================================================================


% ellipsoidal points
center=ell_par(:,1:3)';
a=ell_par(1,4);
b=ell_par(1,5);
c=ell_par(1,6);
alpha=ell_par(1,7);
beta=ell_par(1,8);
gamma=ell_par(1,9);

% scale matrix
A=diag([1/a^2,1/b^2,1/c^2]);

% rotation matrix
invRx=[1 0 0;
    0 cos(-alpha) sin(-alpha);
    0 -sin(-alpha) cos(-alpha);
];

invRy=[cos(-beta) 0 sin(-beta);
    0 1 0;
    -sin(-beta) 0 cos(-beta);
];

invRz=[cos(-gamma) sin(-gamma) 0;
    -sin(-gamma) cos(-gamma) 0;
    0 0 1;
];

R=invRz*invRy*invRx;
M=R'*A*R;
[~, S, V] = svd(M);%LL'=M;
L= real(V' * diag(1./sqrt(diag(S))) * V);



Dimension=3;
NumSamples=num;


% Obtain random samples evenly distributed on the surface of the unit hypersphere
Samples=randn(Dimension,NumSamples);
SampleNorms=sqrt(sum(Samples.^2,1));
Samples=Samples./repmat(SampleNorms,[Dimension 1]); 


% Add some noise
Samples=Samples+0.05*randn(size(Samples));% noise: 0.01-0.05-0.1-0.15-0.2-0.25

% Transform the data into the desired ellipsoid
Samples=L*Samples+repmat(center,[1 NumSamples]); 


Outliers=[];
NumOut=0.3*num;
if type
    Ox=rand(1,NumOut)*200-100;
    Oy=rand(1,NumOut)*200-100;
    Oz=rand(1,NumOut)*150-50;
    Outliers=[Ox;Oy;Oz];
end

Samples1.sample=Samples';
Samples1.outlier=Outliers';

