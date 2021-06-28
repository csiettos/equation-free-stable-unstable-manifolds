% This code finds the stable manifold of the toy disrete-time model
% described in Siettos, C. and Russo, L., 2021, A numerical method for the approximation of stable
% and unstable manifolds of microscopic simulators, Numerical Algorithms,
% DOI: 10.1007/s11075-021-01155-0
c0=zeros(8,1);
x01=0.0;
x02=0.0;
tol=1E-04;
error1=1;
eps=0.01;
ksteps=3;
Jac(1:length(c0),1:length(c0))=0;
while error1>tol
fc0=runmaptoymodel(c0,x01,x02,ksteps);
for i=1:length(c0)
    c=c0;
    c(i)=c0(i)+eps;
    fcp=runmaptoymodel(c,x01,x02,ksteps);
    c(i)=c0(i)-eps;
    fcm=runmaptoymodel(c,x01,x02,ksteps);
    Jac(1:length(c0),i)=(fcp-fcm)/(2*eps);
end
Jacr=eye(length(c))-Jac;
gc0=c0-fc0;
dc=-Jacr\gc0;
c0=c0+dc;
error1=norm(dc);
[error1 c0']
end


    


