function [Transform,iter,spend]=ellipsoid_fit_EM(X, Y,outlierness,normal)

opt.outliers=0.5;   % You can use "outlierness" to set this parameter, but 0.5 works well in most cases. 
opt.max_it=200;     % max number of iterations
opt.tol=1e-8;       % tolerance
opt.normal=normal;

% Convert to double type, save Y
X=double(X);  
Y=double(Y); 


% Volume of the bounding box 
MaxX=max(X);
MinX=min(X);
L=MaxX-MinX;
V=L(1)*L(2)*L(3);


opt.v=V;
normal=opt.normal;

% Start fitting
[R, t,iter,spend]=parameter_solving(X,Y, opt.tol, opt.outliers,opt.v); 
Transform.R=R; 
Transform.t=t;
Transform.s=normal.xscale/normal.yscale;
Transform.t=normal.xscale*t+normal.xd'-Transform.s*(Transform.R*normal.yd');
Transform.R=Transform.s*Transform.R; 

