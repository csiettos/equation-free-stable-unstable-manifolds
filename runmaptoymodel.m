function c=runmaptoymodel(c0,x1n,x2n,ksteps)
dx=0.1;
ksave=1;
      for x01=[x1n-dx x1n+dx]
          for x02=[x2n-dx x2n+dx]
              A0=[x02 x02^2 x01 x01*x02 x01*x02^2 x01^2 x01^2*x02 x01^2*x02^2];
              x03=A0*c0;
              clear x1 x2 x3
              x1(1)=x01;
              x2(1)=x02;
              x3(1)=x03;
              for k=1:ksteps
                  x1(k+1)=-0.5*x1(k);
                  x2(k+1)=-0.5*x2(k)+x1(k)^2;
                  x3(k+1)=2*x3(k)+x2(k)^2;
              end
              x1s(ksave:ksave+ksteps,1)=x1';
              x2s(ksave:ksave+ksteps,1)=x2';
              x3s(ksave:ksave+ksteps,1)=x3';
              ksave=ksave+ksteps+1;
              x1(2:3)=[];x2(2:3)=[];x3(2:3)=[];
              y1=x2s;y2=x2s.^2;y3=x1s;y4=x1s.*x2s;y5=x1s.*x2s.^2;y6=x1s.^2;y7=(x1s.^2).*x2s;y8=(x1s.^2).*x2s.^2;
              A=[y1 y2 y3 y4 y5 y6 y7 y8];
              b=x3s;
          end
      end
      %[U,S,V]=svd(A);
      c=inv(A'*A)*A'*b;