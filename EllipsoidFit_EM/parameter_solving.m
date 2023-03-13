function [B, t,iter_num,iter_time]=parameter_solving(X,Y,tol,outliers,vol)
%===============================================
%
%Inputs:
%-----------------------------------------------
%X                 Nx3 array   the fitted data points
%Y                 Mx3 array   the generated spherical data
%tol               const       tolerance 1e-8 in default
%outliers          outlier weight
%vol               the volume of the bounding box X


%Outputs
%------------------------------------------------
%B                  3x3 the solved affine matrix
%t                  3x1 the translation vector
%=================================================
[N, ~]=size(X);
[M, D]=size(Y);

% Initialization of the variance
sigma2=(M*trace(X'*X)+N*trace(Y'*Y)-2*sum(X)*sum(Y)')/(M*N*D);
T=Y;

% Optimization
iter=0; 
ntol=tol+10; 
L=1;

% Start iteration by the EM
s1=clock;
while ntol >tol
    L_old=L;

   [P1,Pt1,PX,PM,L]=precompute(X,T,sigma2,outliers,vol);
    
    ntol=abs((L-L_old)/L);
  
    % Intermediate parameter 

    Np=sum(P1);
    mu_x=X'*Pt1/Np;
    mu_y=Y'*P1/Np;

    

    % Solve for parameters

    B1=PX'*Y-Np*(mu_x*mu_y');
    B2=(Y.*repmat(P1,1,D))'*Y-Np*(mu_y*mu_y');
    B=B1/B2; 

    t=mu_x-B*mu_y;
    
    sigma2=abs(sum(sum(X.^2.*repmat(Pt1,1,D)))- Np*(mu_x'*mu_x) -trace(B1*B'))/(Np*D); 
     
    % Update weight
    temp=sum(PM);
    outliers=temp/(temp+Np);
 
    
    iter=iter+1;
    
    %Update centroids positioins
    T=Y*B'+repmat(t',[M 1]);
  
end
e1=clock;
iter_num(1)=iter;
iter_time(1)=etime(e1,s1);
end
