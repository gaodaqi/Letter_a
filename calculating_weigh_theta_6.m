function [w3,theta3,minnumerror3,n]=calculating_weigh_theta_6(a22,dim,d22,p0,q0,p2,q2,fid);

n0=1;
nmax=150;
p=p2;
q=q2;
a220=a22;
n1=0;
ir=q0/p0;
while (n0<=(log(p2+q2)/log(2)))

       for n=1:nmax,

          if p==1&q==1
             minnumerror31=zeros(1,3);
             w31=a22(1,2:dim+1)-a22(2,2:dim+1);
             w31=w31/sqrt(w31*w31');
             mu=a22(1,2:dim+1)*w31';
             mu(2)=a22(2,2:dim+1)*w31';
             theta31=-(mu(1)+mu(2))/2;
             n=1;
             numerror(n,1)=0;
             numerror1=0;
             numerror2=0;             
             index0=-1;
             break;
          end

          if n==1
             atrans=a22'*a22;
             arank=rank(atrans);
             adet=det(atrans);
             index0=1;

             if abs(adet)<1.0e-10|arank<dim+1
                atrans=atrans+0.00005*eye(dim+1);
                index0=-1;
%                fprintf('det(%d %d)=%g, rank(%d %d)=%g, Matrix(%d,%d) is close to singular\n',r,t,adet,r,t,arank, r,t);
             end
             invmat=inv(atrans)*a22';
          end

          w=invmat*d22;
          w=w/sqrt(w(2:dim+1)'*w(2:dim+1));
%          theta(n)=w(1)-(p-q)/(2*(p+q));

          project=a22(:,2:dim+1)*w(2:dim+1);
          mu(1)=mean(project(1:p));
          mu(2)=mean(project(p+1:p+q));

%//////////////////////////////////////////////////

          p01=0;
          re1=[];
          for i=1:p,
              if (project(i)<mu(1))&(project(i)>mu(2))
                  p01=p01+1;
                  re1(p01)=project(i);
              end
          end

          q01=0;
          re2=[];
          for i=p+1:p+q,
              if (project(i)<mu(1))&(project(i)>mu(2))
                  q01=q01+1;
                  re2(q01)=project(i);
              end
          end

          if p01==0
             maverage(1)=mu(1);
          else
             maverage(1)=mean(re1);
          end 
          if q01==0
             maverage(2)=mu(2);
          else
             maverage(2)=mean(re2);
          end

%//////////////////////////////////////////////////

          theta(n)=-(maverage(1)+maverage(2))/2;

% Reassigned the desired outputs according to theta

          numerror1=0;
          numerror2=0;
          project00=a220(:,2:dim+1)*w(2:dim+1);
          result00=project00+theta(n);

          result0=project+theta(n);

          for i=1:p2,
              if result00(i)<0
                 numerror1=numerror1+1;
              end
          end
          for i=p2+1:p2+q2,
              if result00(i)>0
                 numerror2=numerror2+1;
              end
          end

%%////////////////////////////////////////////////////////////////////////////
        beta=1.00;
        alpha=sqrt(ir)/50;
        pe=numerror1/p0;
        qe=numerror2/q0;
        if ir<=19
           numerror(n,1)=(p0+q0)*(numerror1+numerror2)/(p0+q0);
        elseif ir<=49
           numerror(n,1)=(p0+q0)*(pe+qe)/2;
        elseif ir<=99.00
           if numerror1>0&numerror2>0
              numerror(n,1)=(p0+q0)*(1-sqrt((1-pe)*(1-qe)));
           elseif numerror1>0&numerror2==0
              numerror(n,1)=(p0+q0)*(1-sqrt(0.0001*(1-pe)));
           elseif numerror1==0&numerror2>0
              numerror(n,1)=(p0+q0)*(1-sqrt(0.0001*(1-qe)));
           else
              numerror(n,1)=0;
           end
        else
           if numerror1>0&numerror2>0
              numerror(n,1)=(p0+q0)*((pe+beta*qe-2*beta*pe*qe)/(2-pe-beta*qe)+alpha*pe)/(1+alpha);
           elseif numerror1>0&numerror2==0
              numerror(n,1)=(p0+q0)*((pe+beta*0.0001-2*beta*pe*0.0001)/(2-pe-beta*0.0001)+alpha*pe)/(1+alpha);
           elseif numerror1==0&numerror2>0
              numerror(n,1)=(p0+q0)*((0.0001+beta*qe-2*beta*0.0001*qe)/(2-0.0001-beta*qe)+alpha*pe)/(1+alpha);
           else
              numerror(n,1)=0;
           end
        end
%        if numtr0(1)>0&numtr0(2)>0
%           numerror(n,1)=(p0+q0)*(1/(1+alpha))*((pe+beta*qe-2*beta*pe*qe)/(2-pe-beta*qe)+alpha*pe);
%        elseif numtr0(1)>0&numtr0(2)==0
%           numerror(n,1)=(p0+q0)*(1/(1+alpha))*((pe+beta*0.0001-2*beta*pe*0.0001)/(2-pe-beta*0.0001)+alpha*pe);
%        elseif numtr0(1)==0&numtr0(2)>0
%           numerror(n,1)=(p0+q0)*(1/(1+alpha))*((0.0001+beta*qe-2*beta*0.0001*qe)/(2-0.0001-beta*qe)+alpha*pe);
%        else
%           numerror(n,1)=0;
%        end
%%/////////////////////////////////////////////////////////////////////////////////////////////
          numerror(n,2)=numerror1;
          numerror(n,3)=numerror2;

          if n==1|rem(n,2)==0
             fprintf('numerror11(%d %d)=%g %g %g %g\n',n0,n,numerror(n,1),numerror1,numerror2,index0);
          end
          if n==1|rem(n,2)==0
             fprintf(fid,'numerror11(%d %d)=%g %g %g %g\n',n0,n,numerror(n,1),numerror1,numerror2,index0);
          end

          if numerror(n,1)==0
             minnumerror33=zeros(1,3);
             w33=w(2:dim+1);
             theta33=theta(n);
             minnumerror31=zeros(1,3);
             w31=w(2:dim+1);
             theta31=theta(n);
             break;
          end

          mmin=1000000000;

          for i=1:p,
              if result0(i)>0
                 if result0(i)<mmin
                    mmin=result0(i);
                 end
              end
          end

          for i=1:p,
              d22(i,1)=result0(i);
              if result0(i)<0
                 d22(i,1)=mmin/2;
              end
          end

          mmax=-10000000;

          for i=p+1:p+q,
              if result0(i)<0
                 if mmax<result0(i)
                    mmax=result0(i);
                 end
              end
          end

          for i=p+1:p+q,
              d22(i,1)=result0(i);
              if result0(i)>0
                 d22(i,1)=mmax/2;
              end
          end

% list 5 continuous weights, thresholds and errors:
              
          if n==1
             w0=zeros(5,dim);
             theta0=zeros(5);
             numerror0=zeros(5,3);
          end

          if n<=5
             w0(n,:)=w(2:dim+1);
             theta0(n)=theta(n);
             numerror0(n,:)=numerror(n,:);
          else
             w0(1:4,:)=w0(2:5,:);
             w0(5,:)=w(2:dim+1);
             theta0(1:4)=theta0(2:5);
             theta0(5)=theta(n);
             numerror0(1:4,:)=numerror0(2:5,:);
             numerror0(5,:)=numerror(n,:);
          end

          if (n>=5)&(n<nmax)&(numerror0(1,1)<min(numerror0(2:5,1)))
              minnumerror33=numerror0(1,:);
              w33=w0(1,:);
              theta33=theta0(1);
              break;
          end

          if (n==nmax)
              numerror01=[numerror0(5,1),numerror0(4,1),numerror0(3,1),numerror0(2,1),numerror0(1,1)];
             [numerror0_0,i00]=min(numerror01);
              minnumerror33=numerror0(5-i00+1,:);
              w33=w0(5-i00+1,:);
              theta33=theta0(5-i00+1);
             break;
          end
       end
       if n>1|rem(n,2)==1
          fprintf('numerror11(%d %d)=%g %g %g %d\n',n0,n,numerror(n,1),numerror1,numerror2,index0);
          fprintf(fid,'numerror11(%d %d)=%g %g %g %d\n',n0,n,numerror(n,1),numerror1,numerror2,index0);
       end

       n1=n1+n;
       if numerror(n,1)==0
          break;
       end

       if n0==1
          w31=w33;
          theta31=theta33;
          minnumerror31=minnumerror33;
       end

       if (n0>1)&(minnumerror33(1,1)>=minnumerror31(1,1))
           break;
       end

       if (n0>1)&(minnumerror33(1,1)<minnumerror31(1,1))
           w31=w33;
           theta31=theta33;
           minnumerror31=minnumerror33;
       end

       w12=w33;
       theta12=theta33;
       a1=a22;
       d=d22;

       if ((p/dim>=8)|(q/dim>=8))&(minnumerror31(1,1)>0)
           [a2,p1,q1,d1]=sample_decomposition_5(a1,w12,dim,d,p,q,p2,q2,fid);
           if (p1==p)&(q1==q);
               break;
           end
           if p1<=1|q1<=1
              break;
           end
           a22=a2;
           d22=d1;
           p=p1;
           q=q1;
           n=1;
%           fprintf('the second-time decomposition is done.\n');
       end
       n0=n0+1;
       if (p/dim<8)&(q/dim<8)&(n0<floor(log(p2+q2)/log(2)))
          n0=floor(log(p2+q2)/log(2));
       end
end
n=n1;
w3=w31;
theta3=theta31;
minnumerror3=minnumerror31;
fprintf('(n, minnumerror3)=(%d, %g, %g, %g)\n',n,minnumerror3);
fprintf(fid,'(n, minnumerror3)=(%d, %g, %g, %g)\n',n,minnumerror3);