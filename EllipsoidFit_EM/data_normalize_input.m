function  [X,pre_normal] =data_normalize_input(x)

n=size(x,1);

pre_normal.xd=mean(x);

x=x-repmat(pre_normal.xd,n,1);

pre_normal.xscale=sqrt(sum(sum(x.^2,2))/n);

X=x/pre_normal.xscale;



