function [w1,theta1,minnumerror,n]=calculating_weigh_theta_5(a1,a0,p0,q0,dim,p,q,d,k0,error_numb,fid);
          
nmax=150;
ir=q0/p0;
for n=1:nmax,

    if p==1&q==1
       minnumerror=zeros(1,3);
       w1=a1(1,2:dim+1)-a1(2,2:dim+1);
       w1=w1/sqrt(w1*w1');
       mu=a1(1,2:dim+1)*w1';
       mu(2)=a1(2,2:dim+1)*w1';
       theta1=-(mu(1)+mu(2))/2;
       index0=0;
       break;
    end

    if n==1
       atrans=a1'*a1;
       arank=rank(atrans);
       adet=det(atrans);
       index0=1;

       if abs(adet)<1.0e-10|arank<dim+1
          atrans=atrans+0.00005*eye(dim+1);
          index0=-1;
%          fprintf('det(%d %d)=%g, rank(%d %d)=%g, Matrix(%d,%d) is close to singular\n',r,t,adet,r,t,arank, r,t);
       end
          invmat=inv(atrans)*a1';
    end

    w=invmat*d;
    w=w/sqrt(w(2:dim+1)'*w(2:dim+1));

%    theta(n)=w(1)-(p-q)/(2*(p+q));

    project=a1(:,2:dim+1)*w(2:dim+1);
    mu(1)=mean(project(1:p));
    mu(2)=mean(project(p+1:p+q));

%//////////////////////////////////////////////////
    p1=0;
    re1=[];
    for i=1:p,
        if (project(i)<mu(1))&(project(i)>mu(2))
            p1=p1+1;
            re1(p1)=project(i);
        end
    end

    q1=0;
    re2=[];
    for i=p+1:p+q,
        if (project(i)<mu(1))&(project(i)>mu(2))
            q1=q1+1;
            re2(q1)=project(i);
        end
    end

    if p1==0
       maverage(1)=mu(1);
    else
       maverage(1)=mean(re1);
    end 
    if q1==0
       maverage(2)=mu(2);
    else
       maverage(2)=mean(re2);
    end
%//////////////////////////////////////////////////
    theta(n)=-(maverage(1)+maverage(2))/2;

% Reassigned the desired outputs according to theta

    numerror1=0;
    numerror2=0;
    project0=a0(:,2:dim+1)*w(2:dim+1);
    result0=project0+theta(n);

    for i=1:p0,
        if result0(i)<0
           numerror1=numerror1+1;
        end
    end
    for i=p0+1:p0+q0,
        if result0(i)>0
           numerror2=numerror2+1;
        end
    end
%%///////////////////////////////////////////////////////////////////////////////
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
          numerror(n,1)=(p0+q0)*(1/(1+alpha))*((pe+beta*qe-2*beta*pe*qe)/(2-pe-beta*qe)+alpha*pe);
       elseif numerror1>0&numerror2==0
          numerror(n,1)=(p0+q0)*(1/(1+alpha))*((pe+beta*0.0001-2*beta*pe*0.0001)/(2-pe-beta*0.0001)+alpha*pe);
       elseif numerror1==0&numerror2>0
          numerror(n,1)=(p0+q0)*(1/(1+alpha))*((0.0001+beta*qe-2*beta*0.0001*qe)/(2-0.0001-beta*qe)+alpha*pe);
       else
          numerror(n,1)=0;
       end
    end
%    if numerror1>0&numerror2>0
%       numerror(n,1)=(p0+q0)*(1/(1+alpha))*((pe+beta*qe-2*beta*pe*qe)/(2-pe-beta*qe)+alpha*pe);
%    elseif numerror1>0&numerror2==0
%       numerror(n,1)=(p0+q0)*(1/(1+alpha))*((pe+beta*0.0001-2*beta*pe*0.0001)/(2-pe-beta*0.0001)+alpha*pe);
%    elseif numerror1==0&numerror2>0
%       numerror(n,1)=(p0+q0)*(1/(1+alpha))*((0.0001+beta*qe-2*beta*0.0001*qe)/(2-0.0001-beta*qe)+alpha*pe);
%    else
%       numerror(n,1)=0;
%    end
%%///////////////////////////////////////////////////////////////////////////////
    numerror(n,2)=numerror1;
    numerror(n,3)=numerror2;

    if n==1|rem(n,1)==0
       fprintf('numerror(%d)=%g %d %d %d\n',n,floor(numerror(n,1)),numerror1,numerror2,index0);
    end

    if n==1|rem(n,1)==0
       fprintf(fid,'numerror(%d)=%g %d %d %d\n',n,floor(numerror(n,1)),numerror1,numerror2,index0);
    end

    if numerror(n,1)==0
       minnumerror=zeros(1,3);
       w1=w(2:dim+1);
       theta1=theta(n);
       break;
    end

    if n==1&numerror(n,1)>2*error_numb
       minnumerror=numerror(n,:);
       w1=w(2:dim+1);
       theta1=theta(n);
       break;
    end

    result1=project+theta(n);
    mmin=10000000;

    for i=1:p0,
        if result0(i)>0
           if result0(i)<mmin
              mmin=result0(i);
           end
        end
    end

%    if numerror1==p0
%       mmin=-((3*mu(1)+mu(2))/4+theta(n));
%    end

    for i=1:p,
        d(i,1)=result1(i);
        if result1(i)<0
           d(i,1)=mmin/2;
        end
    end

    mmax=-10000000;

    for i=p0+1:p0+q0,
        if result0(i)<0
           if mmax<result0(i)
              mmax=result0(i);
           end
        end
    end

    for i=p+1:p+q,
        d(i,1)=result1(i);
        if result1(i)>0
           d(i,1)=mmax/2;
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

%if ir<=7
    if (n>=5)&(n<nmax)&(numerror0(1,1)<min(numerror0(2:5,1)))
        minnumerror=numerror0(1,:);
        w1=w0(1,:);
        theta1=theta0(1);
        break;
    end
%else
%    if (n>=5)&(n<nmax)&(numerror0(1,1)<min(numerror0(2:5,1)))&(numerror0(1,2)<min(numerror0(2:5,2)))
%        minnumerror=numerror0(1,:);
%        w1=w0(1,:);
%        theta1=theta0(1);
%        break;
%    end
%end

    if (n==nmax)
        numerror01=[numerror0(5,1),numerror0(4,1),numerror0(3,1),numerror0(2,1),numerror0(1,1)];
       [numerror0_0,i00]=min(numerror01);
        minnumerror=numerror0(5-i00+1,:);
        w1=w0(5-i00+1,:);
        theta1=theta0(5-i00+1);
        break;
    end
end
%fprintf('numerror(%d)=%d %d %d %d\n',n,numerror(n,1),numerror1,numerror2,index0);
%fprintf(fid,'numerror(%d)=%d %d %d %d\n',n,numerror(n,1),numerror1,numerror2,index0);