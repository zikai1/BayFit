%===========================
% Robust Ellipsoid-specific Fitting via expectation-maximization 
% Update 3.13.2023
%============================

clc;
close all;
clear;


%% Test example ellipsoid
ellParFit=[5 7 18 35 65 90 80/180*pi 50/180*pi 0/180*pi];


%% Generate example ellipsoid
ptFit=generate_ellipsoidal_data(ellParFit,1,200);

%% Compute outlierness
ptFit=[ptFit.sample;ptFit.outlier];
[rdos_score,X,X_normal]=outlier_det(ptFit);
inlier=ptFit(rdos_score<=2,:);
inlier_num=size(inlier,1);
init_center=mean(inlier,1);% mass center
outlierness=1-inlier_num/size(ptFit,1);

%% Initialization the ellipsoid parameter
ellParInit=[init_center 1 1 1 0 0 0];% init_center
ptInit=drawEllipsoid(ellParInit,1,sqrt(inlier_num));   

%% Normalization 
[Y,Y_normal]=data_normalize_input(ptInit);

%% Start fitting
normal.xd=X_normal.xd;
normal.yd=Y_normal.xd;
normal.xscale=X_normal.xscale;
normal.yscale=Y_normal.xscale;
[transform,iter,spend]=ellipsoid_fit_EM(X,Y,outlierness,normal);

%% Get final ellipsoid   
[FitEllipsoid]=ellipsoid_par(init_center,transform.R,transform.t);
    
%% Plot result   
figure
subplot(1,2,1);
plot3(ptFit(:,1),ptFit(:,2),ptFit(:,3),'r.');
title("Input data");

subplot(1,2,2);
hold on;
plot3(ptFit(:,1),ptFit(:,2),ptFit(:,3),'r.');
Fit=drawEllipsoid(FitEllipsoid,2,40);
title("The fitted ellipsoid");