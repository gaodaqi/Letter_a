clear all;

t0=clock;

load D:\highimbalance\letreca\data_5rand\rand_train1.txt;
load D:\highimbalance\letreca\ntr.txt;
dim=size(rand_train1,2)-2;
%dim=8;
sign0=[];
k0=1;
c=size(ntr,1);
fid=fopen('D:\highimbalance\letreca\inv_dc1\result\basicdata_rand_train1.txt','w');
total_time=zeros(1,c*(c-1)/2);

for r=1:c-1,
    for t=r+1:c,

        p0=ntr(r,2)-ntr(r,1)+1;
        q0=ntr(t,2)-ntr(t,1)+1;
        ir=q0/p0;
        d00=ones(p0,1);
        d00(p0+1:p0+q0,1)=-ones(q0,1);
        d0=d00;        

        d=d0;
        a0=[];
        a0(:,2:dim+1)=rand_train1(ntr(r,1):ntr(r,2),1:dim);

        if p0==1
           mu01=a0(:,2:dim+1);
        else
           mu01=mean(a0(:,2:dim+1));
        end
        a0(p0+1:p0+q0,2:dim+1)=rand_train1(ntr(t,1):ntr(t,2),1:dim);        
        if q0==1
           mu02=a0(p0+1:p0+1,2:dim+1);
        else
           mu02=mean(a0(p0+1:p0+q0,2:dim+1));
        end
        a0(:,1)=ones(p0+q0,1);

        a1=a0;

        p=p0;
        q=q0;

        n0=1;
        n000(k0)=0;
        sign0(k0)=1;
        error_numb=100000000;

        while (n0<=(log(p0+q0)/log(2)))
              [w1,theta1,minnumerror,n]=calculating_weigh_theta_5(a1,a0,p0,q0,dim,p,q,d,k0,error_numb,fid);
               minnumerror1=floor(10*minnumerror(1));
               minnumerror1_1=floor(minnumerror(2));
               n000(k0)=n000(k0)+n;
               if n0==1
                  w2(k0,:)=w1';
                  theta2(k0)=theta1;
                  minnumerror2(k0)=minnumerror1;
                  minnumerror2_1(k0)=minnumerror1_1;
                  error_numb=minnumerror1;
               end
               if n0==1&minnumerror1==0
                   fprintf('n0=%d;k0=%d;minnumerror1=%g;minnumerror2(k0)=%g\n',n0,k0,minnumerror1,minnumerror2(k0));
                   fprintf(fid,'n0=%d;k0=%d;minnumerror1=%g;minnumerror2(k0)=%g\n',n0,k0,minnumerror1,minnumerror2(k0));
                   break;
               end
               if (n0>1)&(minnumerror1==minnumerror2(k0))&(floor(sum(w1==w2(k0,:))/dim))&(theta1==theta2(k0))
                  break;
               end
               if (n0>1)&(minnumerror1>1.5*minnumerror2(k0))
                  break;
               end
               if (n0>1)&(minnumerror1<minnumerror2(k0))&(((ir>7)&minnumerror1_1<=minnumerror2_1(k0))|(ir<=7))
                   w2(k0,:)=w1';
                   theta2(k0)=theta1;
                   minnumerror2(k0)=minnumerror1;
                   minnumerror2_1(k0)=minnumerror1_1;
                   error_numb=minnumerror1;
               end
               if (n0>1)&(minnumerror1==0|minnumerror1>1.5*minnumerror2(k0))
                   fprintf('n0=%d;k0=%d;minnumerror1=%g;minnumerror2(k0)=%g\n',n0,k0,minnumerror1,minnumerror2(k0));
                   fprintf(fid,'n0=%d;k0=%d;minnumerror1=%g;minnumerror2(k0)=%g\n',n0,k0,minnumerror1,minnumerror2(k0));
                   break;
               end
               if (p/dim<8)&(q/dim<8)
                   fprintf('n0=%d;k0=%d;minnumerror1=%g;minnumerror2(k0)=%g\n',n0,k0,minnumerror1,minnumerror2(k0));
                   fprintf(fid,'n0=%d;k0=%d;minnumerror1=%g;minnumerror2(k0)=%g\n',n0,k0,minnumerror1,minnumerror2(k0));
                   break;
               end
               w12=w2(k0,:);
               theta12=theta2(k0);

               if ((p/dim>=8)|(q/dim>=8))
                   [a2,p1,q1,d1]=sample_decomposition_5(a1,w12,dim,d,p,q,p0,q0,fid);
                   if (p1==p)&(q1==q);
                       break;
                   end
                   a1=a2;
                   d=d1;
                   p=p1;
                   q=q1;
               end
               if p<=1|q<=1
                  break;
               end
               n0=n0+1;
        end

        fprintf('num(%d,%d)=%d\n',r,t,n000(k0));
        fprintf(fid,'num(%d,%d)=%d\n',r,t,n000(k0));
        w12=w2(k0,:);
        theta12=theta2(k0);

        numerror1=0;
        numerror2=0;
        project0=a0(:,2:dim+1)*w12';
        result0=project0+theta12;

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
        fprintf('minnumerror02(%d,%d)=%d  %d  %d\n',r,t,numerror1,numerror2,numerror1+numerror2);
        fprintf(fid,'minnumerror02(%d,%d)=%d  %d  %d\n',r,t,numerror1,numerror2,numerror1+numerror2);

        numtr0=[0 0];
        if sign0(k0)==1
           result02=a0(:,2:dim+1)*w12'+ones(p0+q0,1)*theta12;
           re_mu(1)=mu01*w12'+theta12;
           re_mu(2)=mu02*w12'+theta12;
           numtr1=zeros(4,2);
           for i=1:p0,
               if result02(i)<0&result02(i)>=re_mu(2)
                  numtr0(1)=numtr0(1)+1;
                  numtr1(2,1)=numtr1(2,1)+1;
               end
               if result02(i)<0&result02(i)<re_mu(2)
                  numtr0(1)=numtr0(1)+1;
                  numtr1(4,1)=numtr1(4,1)+1;
               end
           end
           for i=p0+1:p0+q0,
               if result02(i)>=0&result02(i)<=re_mu(1)
                  numtr0(2)=numtr0(2)+1;
                  numtr1(1,2)=numtr1(1,2)+1;
               end
               if result02(i)>=0&result02(i)>re_mu(1)
                  numtr0(2)=numtr0(2)+1;
                  numtr1(3,2)=numtr1(3,2)+1;
               end
           end
           numtr0
           numtr1'
           sum(sum(numtr1))
           minnumerror2(k0)=sum(numtr0);
           fprintf(fid,'%g  %g\n',numtr0');
           fprintf(fid,'%g  %g\n',numtr1');
           fprintf(fid,'%g\n',sum(sum(numtr1)));
           minnumerror01=[];
           minnumerror01=sum(numtr1');
        end
        t1=etime(clock,t0)
        fprintf(fid,'runtime(%d,%d)=%g\n',r,t,t1);
        aver_index=[];
        numtr_0=[];
        tacc_index=[];
        numtr_0(1)=numtr0(1);
        tacc_index(1)=sum(numtr0);
%%////////////////////////////////////////////////////////////////////////////
        beta=1.00;
        alpha=sqrt(ir)/50;
        pe=numtr0(1)/p0;
        qe=numtr0(2)/q0;
        if ir<=19
           aver_index(1)=(numtr0(1)+numtr0(2))/(p0+q0);
        elseif ir<=49
           aver_index(1)=(pe+qe)/2;
        elseif ir<=99.00
           if numtr0(1)>0&numtr0(2)>0
              aver_index(1)=(1-sqrt((1-pe)*(1-qe)));
           elseif numtr0(1)>0&numtr0(2)==0
              aver_index(1)=(1-sqrt(0.0001*(1-pe)));
           elseif numtr0(1)==0&numtr0(2)>0
              aver_index(1)=(1-sqrt(0.0001*(1-qe)));
           else
              aver_index(1)=0;
           end
        else
           if numtr0(1)>0&numtr0(2)>0
              aver_index(1)=(1/(1+alpha))*((pe+beta*qe-2*beta*pe*qe)/(2-pe-beta*qe)+alpha*pe);
           elseif numtr0(1)>0&numtr0(2)==0
              aver_index(1)=(1/(1+alpha))*((pe+beta*0.0001-2*beta*pe*0.0001)/(2-pe-beta*0.0001)+alpha*pe);
           elseif numtr0(1)==0&numtr0(2)>0
              aver_index(1)=(1/(1+alpha))*((0.0001+beta*qe-2*beta*0.0001*qe)/(2-0.0001-beta*qe)+alpha*pe);
           else
              aver_index(1)=0;
           end
        end
%%/////////////////////////////////////////////////////////////////////////////////////////////
        if aver_index(1)>0.01
           [a3,d3,n_tr3]=sample_decomposition_6(a0,w12,theta12,dim,p0,q0);
            n_tr3

            fprintf(fid,'\n n_tr3=\n\n',n_tr3');
            fprintf(fid,'%g  %g  %g  %g\n',n_tr3');
            fprintf(fid,'\n');
            minnumerror3_1=zeros(1,4);
            minnumerror3_2=zeros(4,2);
            for i=1,
                if n_tr3(2*i-1,3)==0&n_tr3(2*i,3)>0
                   w3_1(i,:)=-w2(k0,:);
                   theta3_1(i)=-theta2(k0);
                   n001(i)=0;
                end
                if n_tr3(2*i-1,3)>=0&n_tr3(2*i,3)==0
                   w3_1(i,:)=w2(k0,:);
                   theta3_1(i)=theta2(k0);
                   n001(i)=0;
                end
            end
            for i=2,
                if n_tr3(2*i-1,3)==0&n_tr3(2*i,3)>=0
                   w3_1(i,:)=w2(k0,:);
                   theta3_1(i)=theta2(k0);
                   n001(i)=0;
                end
                if n_tr3(2*i-1,3)>0&n_tr3(2*i,3)==0
                   w3_1(i,:)=-w2(k0,:);
                   theta3_1(i)=-theta2(k0);
                   n001(i)=0;
                end
            end
            for i=3
                if n_tr3(2*i-1,3)==0&n_tr3(2*i,3)>0
                   w3_1(i,:)=-w2(k0,:);
                   theta3_1(i)=mu01*w2(k0,:)';
                   n001(i)=0;
                end
                if n_tr3(2*i-1,3)>=0&n_tr3(2*i,3)==0
                   w3_1(i,:)=w2(k0,:);
                   theta3_1(i)=-mu01*w2(k0,:)';
                   n001(i)=0;
                end
            end
            for i=4
                if n_tr3(2*i-1,3)==0&n_tr3(2*i,3)>=0
                   w3_1(i,:)=w2(k0,:);
                   theta3_1(i)=-mu02*w2(k0,:)';
                   n001(i)=0;
                end
                if n_tr3(2*i-1,3)>0&n_tr3(2*i,3)==0
                   w3_1(i,:)=-w2(k0,:);
                   theta3_1(i)=mu02*w2(k0,:)';
                   n001(i)=0;
                end
            end

            for ii=1:4,
                p2=n_tr3(2*ii-1,3);
                q2=n_tr3(2*ii,3);
                a22=[];
                d22=[];
                if p2>0&q2>0
                   a22=a3(n_tr3(2*ii-1,1):n_tr3(2*ii,2),:);
                   d22=d3(n_tr3(2*ii-1,1):n_tr3(2*ii,2),:);
                  [w3,theta3,minnumerror3,n]=calculating_weigh_theta_6(a22,dim,d22,p0,q0,p2,q2,fid);
%                   minnumerror3_0=floor(10*minnumerror3(1));
                   fprintf('num(%d,%d,%d)=%d\n\n',r,t,ii,n);
                   fprintf(fid,'num(%d,%d,%d)=%d\n\n',r,t,ii,n);
                   w3_1(ii,:)=w3';
                   theta3_1(ii)=theta3;
                   minnumerror3_1(ii)=sum(minnumerror3(2:3));
                   minnumerror3_2(ii,:)=minnumerror3(2:3);
                   n001(ii)=n;
                end
            end
            n000(k0)=n000(k0)+sum(n001);

            s=sum(sign0);
            w4(s,:)=w2(k0,:);
            theta4(s)=theta2(k0);

            fprintf('minnumerror3_1(1,2,3,4)=%d %d %d %d, sum=%g\n',minnumerror3_1(1),minnumerror3_1(2),minnumerror3_1(3),minnumerror3_1(4),sum(minnumerror3_1));
            fprintf(fid,'minnumerror3_1(1,2,3,4)=%d %d %d %d, sum=%g\n',minnumerror3_1(1),minnumerror3_1(2),minnumerror3_1(3),minnumerror3_1(4),sum(minnumerror3_1));
               sign0(k0)=7;
               minnumerror4_1=[0 0 0 0];

               if minnumerror3_1(1)<min([n_tr3(2,3),minnumerror01(1)])&minnumerror3_2(1,2)<n_tr3(2,3)
                  w4_1(1,:)=w3_1(1,:);
                  theta4_1(1)=theta3_1(1);
                  minnumerror4_1(1)=minnumerror3_1(1);
               else
                  w4_1(1,:)=w2(k0,:);
                  theta4_1(1)=theta2(k0);
                  minnumerror4_1(1)=n_tr3(2,3);
               end

               if minnumerror3_1(2)<min([n_tr3(3,3),minnumerror01(2)])&minnumerror3_2(2,1)<n_tr3(3,3)
                  w4_1(2,:)=w3_1(2,:);
                  theta4_1(2)=theta3_1(2);
                  minnumerror4_1(2)=minnumerror3_1(2);
               else
                  w4_1(2,:)=w2(k0,:);
                  theta4_1(2)=theta2(k0);
                  minnumerror4_1(2)=n_tr3(3,3);
               end

               if minnumerror3_1(3)<min([n_tr3(6,3),minnumerror01(3)])&minnumerror3_2(3,2)<n_tr3(6,3)
                  w4_1(3,:)=w3_1(3,:);
                  theta4_1(3)=theta3_1(3);
                  minnumerror4_1(3)=minnumerror3_1(3);
               else
                  w4_1(3,:)=w2(k0,:);
                  theta4_1(3)=-mu01*w2(k0,:)';
                  minnumerror4_1(3)=n_tr3(6,3);
               end

               if minnumerror3_1(4)<min([n_tr3(7,3),minnumerror01(4)])&minnumerror3_2(4,1)<n_tr3(7,3)
                  w4_1(4,:)=w3_1(4,:);
                  theta4_1(4)=theta3_1(4);
                  minnumerror4_1(4)=minnumerror3_1(4);
               else
                  w4_1(4,:)=w2(k0,:);
                  theta4_1(4)=-mu02*w2(k0,:)';
                  minnumerror4_1(4)=n_tr3(7,3);
               end

               w4(s+1,:)=w2(k0,:);
               theta4(s+1)=-mu01*w2(k0,:)';
               w4(s+2,:)=w2(k0,:);
               theta4(s+2)=-mu02*w2(k0,:)';
               w4(s+3:s+6,:)=w4_1;
               theta4(s+3:s+6)=theta4_1;

               minnumerror4=minnumerror4_1;
               minnumerror2(k0)=sum(minnumerror4);
               fprintf('minnumerror2(%d,%d)=minnumerror2(%d)=%g\n',r,t,k0,minnumerror2(k0));
               fprintf(fid,'minnumerror2(%d,%d)=minnumerror2(%d)=%g\n',r,t,k0,minnumerror2(k0));
               fprintf('minnumerror4(1, 2, 3, 4)=%g %g %g %g, sum=%g\n',minnumerror4(1),minnumerror4(2),minnumerror4(3),minnumerror4(4),sum(minnumerror4));
               fprintf(fid,'minnumerror4(1, 2, 3, 4)=%g %g %g %g, sum=%g\n',minnumerror4(1),minnumerror4(2),minnumerror4(3),minnumerror4(4),sum(minnumerror4));

               w44=[];
               theta44=[];
               w44=w4_1;
               theta44=theta4_1;
%%/////////////////////////////////////////////////////////////////////////////////////////////
               kk=sum(sign0(1:k0));
               numtr0=[0 0];
               if sign0(k0)==7
                  re=a0(:,2:dim+1)*(w4(kk-6:kk,:))'+ones(p0+q0,1)*theta4(kk-6:kk);
                  numtr1=zeros(4,2);
                  for i=1:p0,
                      if (re(i,1)>0)&(re(i,2)<0)&(re(i,4)<0)
                          numtr0(1)=numtr0(1)+1;
                          numtr1(1,1)=numtr1(1,1)+1;
                      end
                      if (re(i,1)<0)&(re(i,3)>0)&(re(i,5)<0)
                          numtr0(1)=numtr0(1)+1;
                          numtr1(2,1)=numtr1(2,1)+1;
                      end
                      if (re(i,1)>0)&(re(i,2)>0)&(re(i,6)<0)
                          numtr0(1)=numtr0(1)+1;
                          numtr1(3,1)=numtr1(3,1)+1;
                      end
                      if (re(i,1)<0)&(re(i,3)<0)&(re(i,7)<0)
                          numtr0(1)=numtr0(1)+1;
                          numtr1(4,1)=numtr1(4,1)+1;
                      end
                  end
                  for i=p0+1:p0+q0,
                      if (re(i,1)>0)&(re(i,2)<0)&(re(i,4)>0)
                          numtr0(2)=numtr0(2)+1;
                          numtr1(1,2)=numtr1(1,2)+1;
                      end
                      if (re(i,1)<0)&(re(i,3)>0)&(re(i,5)>0)
                          numtr0(2)=numtr0(2)+1;
                          numtr1(2,2)=numtr1(2,2)+1;
                      end
                      if (re(i,1)>0)&(re(i,2)>0)&(re(i,6)>0)
                          numtr0(2)=numtr0(2)+1;
                          numtr1(3,2)=numtr1(3,2)+1;
                      end
                      if (re(i,1)<0)&(re(i,3)<0)&(re(i,7)>0)
                          numtr0(2)=numtr0(2)+1;
                          numtr1(4,2)=numtr1(4,2)+1;
                      end
                  end
                  numtr0
                  numtr1'
                  sum(sum(numtr1))
                  minnumerror2(k0)=sum(numtr0);
                  minnumerror4=sum(numtr1');
                  fprintf(fid,'%g  %g\n',numtr0');
                  fprintf(fid,'%g  %g\n',numtr1');
                  fprintf(fid,'%g\n',sum(sum(numtr1)));
                  fprintf('minnumerror2(%d,%d)=minnumerror2(%d)=%g\n',r,t,k0,minnumerror2(k0));
                  fprintf(fid,'minnumerror2(%d,%d)=minnumerror2(%d)=%g\n',r,t,k0,minnumerror2(k0));
                  fprintf('minnumerror4(1, 2, 3, 4)=%g %g %g %g, sum=%g\n',minnumerror4(1),minnumerror4(2),minnumerror4(3),minnumerror4(4),sum(minnumerror4));
                  fprintf(fid,'minnumerror4(1, 2, 3, 4)=%g %g %g %g, sum=%g\n',minnumerror4(1),minnumerror4(2),minnumerror4(3),minnumerror4(4),sum(minnumerror4));
               end
         end
         numtr_0(2)=numtr0(1);
         tacc_index(2)=sum(numtr0);

%%////////////////////////////////////////////////////////////////////////////
        pe=numtr0(1)/p0;
        qe=numtr0(2)/q0;
        if ir<=19
           aver_index(2)=(numtr0(1)+numtr0(2))/(p0+q0);
        elseif ir<=49
           aver_index(2)=(pe+qe)/2;
        elseif ir<=99.00
           if numtr0(1)>0&numtr0(2)>0
              aver_index(2)=(1-sqrt((1-pe)*(1-qe)));
           elseif numtr0(1)>0&numtr0(2)==0
              aver_index(2)=(1-sqrt(0.0001*(1-pe)));
           elseif numtr0(1)==0&numtr0(2)>0
              aver_index(2)=(1-sqrt(0.0001*(1-qe)));
           else
              aver_index(2)=0;
           end
        else
           if numtr0(1)>0&numtr0(2)>0
              aver_index(2)=(1/(1+alpha))*((pe+beta*qe-2*beta*pe*qe)/(2-pe-beta*qe)+alpha*pe);
           elseif numtr0(1)>0&numtr0(2)==0
              aver_index(2)=(1/(1+alpha))*((pe+beta*0.0001-2*beta*pe*0.0001)/(2-pe-beta*0.0001)+alpha*pe);
           elseif numtr0(1)==0&numtr0(2)>0
              aver_index(2)=(1/(1+alpha))*((0.0001+beta*qe-2*beta*0.0001*qe)/(2-0.0001-beta*qe)+alpha*pe);
           else
              aver_index(2)=0;
           end
        end
%%/////////////////////////////////////////////////////////////////////////////////////////////
         t1=etime(clock,t0)
         fprintf(fid,'runtime(%d,%d)=%g\n',r,t,t1);
%/////////////////////////////////////////////////////////////////////////////////
         if (aver_index(2)>0.025)&(((ir>99)&(numtr_0(2)<=numtr_0(1)))|(ir<=99))
         
           [w5,theta5]=calculating_weigh_theta_7(mu01,mu02,w12,theta12,dim,fid);

            a5=[];
            d5=[];
            a5_1=[];
            d5_1=[];
            a5_2=[];
            d5_2=[];
            n_tr5=[];
            n_tr5_1=[];
            n_tr5_2=[];

            n_tr5_1=zeros(8,3);
            n_tr5_2=zeros(8,3);
            if n_tr3(1,3)==0&n_tr3(2,3)>0
               n_tr5_1(1,1:3)=zeros(1,3);
               n_tr5_1(2,1)=1;
               n_tr5_1(2,2)=floor(n_tr3(2,3)/2);
               n_tr5_1(2,3)=n_tr5_1(2,2);
               a5_1=a3(n_tr3(2,1):n_tr3(2,1)-1+floor(n_tr3(2,3)/2),:);
               d5_1(:,1)=-ones(floor(n_tr3(2,3)/2),1);

               n_tr5_2(1,1:3)=zeros(1,3);
               n_tr5_2(2,1)=1;
               n_tr5_2(2,2)=n_tr3(2,3)-floor(n_tr3(2,3)/2);
               n_tr5_2(2,3)=n_tr5_2(2,2);
               a5_2=a3(n_tr3(2,1)+floor(n_tr3(2,3)/2)+1:n_tr3(2,2),:);
               d5_2(:,1)=-ones(n_tr3(2,3)-floor(n_tr3(2,3)/2),1);
            end
            if n_tr3(1,3)>0&n_tr3(2,3)==0
               if n_tr3(1,3)>1
                  n_tr5_1(1,1)=1;
                  n_tr5_1(1,2)=floor(n_tr3(1,3)/2);
                  n_tr5_1(1,3)=n_tr5_1(1,2);
                  a5_1=a3(n_tr3(1,1):floor(n_tr3(1,3)/2),:);
                  d5_1(:,1)=ones(floor(n_tr3(1,3)/2),1);
               else
                  n_tr5_1(1,1:3)=zeros(1,3);
               end
               n_tr5_1(2,1:3)=zeros(1,3);
               n_tr5_2(1,1)=1;
               n_tr5_2(1,2)=n_tr3(1,3)-floor(n_tr3(1,3)/2);
               n_tr5_2(1,3)=n_tr5_2(1,2);
               n_tr5_2(2,1:3)=zeros(1,3);
               a5_2=a3(floor(n_tr3(1,3)/2)+1:n_tr3(1,3),:);
               d5_2(:,1)=ones(n_tr3(1,3)-floor(n_tr3(1,3)/2),1);
            end
            if (n_tr3(1,3)>0&n_tr3(2,3)>0)&minnumerror4(1)==0
               if n_tr3(1,3)>1
                  n_tr5_1(1,1)=1;
                  n_tr5_1(1,2)=floor(n_tr3(1,3)/2);
                  n_tr5_1(1,3)=n_tr5_1(1,2);
                  a5_1=a3(n_tr3(1,1):floor(n_tr3(1,3)/2),:);
                  d5_1(:,1)=ones(floor(n_tr3(1,3)/2),1);
               else
                  n_tr5_1(1,1:3)=zeros(1,3);
               end
               n_tr5_1(2,1)=n_tr5_1(1,2)+1;
               n_tr5_1(2,2)=n_tr5_1(1,2)+floor(n_tr3(2,3)/2);
               n_tr5_1(2,3)=floor(n_tr3(2,3)/2);

               a5_1(n_tr5_1(2,1):n_tr5_1(2,2),:)=a3(n_tr3(2,1):n_tr3(2,1)+floor(n_tr3(2,3)/2)-1,:);
               d5_1(n_tr5_1(2,1):n_tr5_1(2,2),1)=-ones(floor(n_tr3(2,3)/2),1);

               n_tr5_2(1,1)=1;
               n_tr5_2(1,2)=n_tr3(1,3)-floor(n_tr3(1,3)/2);
               n_tr5_2(1,3)=n_tr5_2(1,2);
               n_tr5_2(2,1)=n_tr5_2(1,2)+1;
               n_tr5_2(2,2)=n_tr5_2(1,2)+n_tr3(2,3)-floor(n_tr3(2,3)/2);
               n_tr5_2(2,3)=n_tr3(2,3)-floor(n_tr3(2,3)/2);

               a5_2=a3(floor(n_tr3(1,3)/2)+1:n_tr3(1,3),:);
               a5_2(n_tr5_2(2,1):n_tr5_2(2,2),:)=a3(n_tr3(1,2)+floor(n_tr3(2,3)/2)+1:n_tr3(2,2),:);
               d5_2(:,1)=ones(n_tr3(1,3)-floor(n_tr3(1,3)/2),1);
               d5_2(n_tr5_2(2,1):n_tr5_2(2,2),1)=-ones(n_tr3(2,3)-floor(n_tr3(2,3)/2),1);
            end
            if n_tr3(1,3)>0&n_tr3(2,3)>0&minnumerror4(1)>0
               n_tr4=[];
               n_tr4(1,1)=1;
               n_tr4(1,2)=n_tr3(1,3);
               n_tr4(1,3)=n_tr3(1,3);
               n_tr4(1,4)=1;
               n_tr4(2,1)=n_tr3(1,3)+1;
               n_tr4(2,2)=n_tr3(1,3)+n_tr3(2,3);
               n_tr4(2,3)=n_tr3(2,3);
               n_tr4(2,4)=2;
               a00=a3(n_tr3(1,1):n_tr3(2,2),:);

              [a51,d51,n_tr51]=sample_decomposition_7(a00,n_tr4,w5,theta5,dim);
               n_tr5_1(1:2,1:3)=n_tr51(1:2,1:3);

               if n_tr51(1,3)==0&n_tr51(2,3)>0
                  a5_1=a51(n_tr51(2,1):n_tr51(2,2),:);
                  d5_1=d51(n_tr51(2,1):n_tr51(2,2),:);
               end
               if n_tr51(1,3)>0&n_tr51(2,3)==0
                  a5_1=a51(n_tr51(1,1):n_tr51(1,2),:);
                  d5_1=d51(n_tr51(1,1):n_tr51(1,2),:);
               end
               if n_tr51(1,3)>0&n_tr51(2,3)>0
                  a5_1=a51(n_tr51(1,1):n_tr51(2,2),:);
                  d5_1=d51(n_tr51(1,1):n_tr51(2,2),:);
               end

               if n_tr51(3,3)==0&n_tr51(4,3)>0
                  n_tr5_2(1,1:3)=zeros(1,3);
                  n_tr5_2(2,1)=1;
                  n_tr5_2(2,2)=n_tr51(4,3);
                  n_tr5_2(2,3)=n_tr51(4,3);
                  a5_2=a51(n_tr51(4,1):n_tr51(4,2),:);
                  d5_2=d51(n_tr51(4,1):n_tr51(4,2),:);
               end
               if n_tr51(3,3)>0&n_tr51(4,3)==0
                  n_tr5_2(1,1)=1;
                  n_tr5_2(1,2)=n_tr51(3,3);
                  n_tr5_2(1,3)=n_tr51(3,3);
                  n_tr5_2(2,1:3)=zeros(1,3);
                  a5_2=a51(n_tr51(3,1):n_tr51(3,2),:);
                  d5_2=d51(n_tr51(3,1):n_tr51(3,2),:);
               end
               if n_tr51(3,3)>0&n_tr51(4,3)>0
                  n_tr5_2(1,1)=1;
                  n_tr5_2(1,2)=n_tr51(3,3);
                  n_tr5_2(1:2,3)=n_tr51(3:4,3);
                  n_tr5_2(2,1)=n_tr51(3,3)+1;
                  n_tr5_2(2,2)=n_tr51(3,3)+n_tr51(4,3);
                  a5_2=a51(n_tr51(3,1):n_tr51(4,2),:);
                  d5_2=d51(n_tr51(3,1):n_tr51(4,2),:);
               end
            end

            for ii=2:4
                if n_tr3(2*ii-1,3)==0&n_tr3(2*ii,3)>0
                   n_tr5_1(2*ii-1,1:3)=zeros(1,3);
                   n_tr5_1(2*ii,1)=sum(n_tr5_1(1:2*(ii-1),3))+1;
                   n_tr5_1(2*ii,2)=sum(n_tr5_1(1:2*(ii-1),3))+floor(n_tr3(2*ii,3)/2);
                   n_tr5_1(2*ii,3)=floor(n_tr3(2*ii,3)/2);
                   a5_1(n_tr5_1(2*ii,1):n_tr5_1(2*ii,2),:)=a3(sum(n_tr3(1:2*(ii-1),3))+1:sum(n_tr3(1:2*(ii-1),3))+floor(n_tr3(2*ii,3)/2),:);
                   d5_1(n_tr5_1(2*ii,1):n_tr5_1(2*ii,2),1)=-ones(floor(n_tr3(2*ii,3)/2),1);

                   n_tr5_2(2*ii-1,1:3)=zeros(1,3);
                   n_tr5_2(2*ii,1)=sum(n_tr5_2(1:2*(ii-1),3))+1;
                   n_tr5_2(2*ii,2)=sum(n_tr5_2(1:2*(ii-1),3))+n_tr3(2*ii,3)-floor(n_tr3(2*ii,3)/2);
                   n_tr5_2(2*ii,3)=n_tr3(2*ii,3)-floor(n_tr3(2*ii,3)/2);
                   a5_2(n_tr5_2(2*ii,1):n_tr5_2(2*ii,2),:)=a3(sum(n_tr3(1:2*(ii-1),3))+floor(n_tr3(2*ii,3)/2)+1:sum(n_tr3(1:2*(ii-1),3))+n_tr3(2*ii,3),:);
                   d5_2(n_tr5_2(2*ii,1):n_tr5_2(2*ii,2),1)=-ones(n_tr3(2*ii,3)-floor(n_tr3(2*ii,3)/2),1);
                end
                if n_tr3(2*ii-1,3)>0&n_tr3(2*ii,3)==0
                   n_tr5_1(2*ii-1,1)=sum(n_tr5_1(1:2*(ii-1),3))+1;
                   n_tr5_1(2*ii-1,2)=sum(n_tr5_1(1:2*(ii-1),3))+floor(n_tr3(2*ii-1,3)/2);
                   n_tr5_1(2*ii-1,3)=floor(n_tr3(2*ii-1,3)/2);
                   a5_1(n_tr5_1(2*ii-1,1):n_tr5_1(2*ii-1,2),:)=a3(sum(n_tr3(1:2*(ii-1),3))+1:sum(n_tr3(1:2*(ii-1),3))+floor(n_tr3(2*ii-1,3)/2),:);
                   d5_1(n_tr5_1(2*ii-1,1):n_tr5_1(2*ii-1,2),1)=ones(floor(n_tr3(2*ii-1,3)/2),1);
                   n_tr5_1(2*ii,1:3)=zeros(1,3);

                   n_tr5_2(2*ii-1,1)=sum(n_tr5_2(1:2*(ii-1),3))+1;
                   n_tr5_2(2*ii-1,2)=sum(n_tr5_2(1:2*(ii-1),3))+n_tr3(2*ii-1,3)-floor(n_tr3(2*ii-1,3)/2);
                   n_tr5_2(2*ii-1,3)=n_tr3(2*ii-1,3)-floor(n_tr3(2*ii-1,3)/2);
                   a5_2(n_tr5_2(2*ii-1,1):n_tr5_2(2*ii-1,2),:)=a3(sum(n_tr3(1:2*(ii-1),3))+floor(n_tr3(2*ii-1,3)/2)+1:sum(n_tr3(1:2*(ii-1),3))+n_tr3(2*ii-1,3),:);
                   d5_2(n_tr5_2(2*ii-1,1):n_tr5_2(2*ii-1,2),1)=ones(n_tr3(2*ii-1,3)-floor(n_tr3(2*ii-1,3)/2),1);
                   n_tr5_2(2*ii,1:3)=zeros(1,3);
                end

                if n_tr3(2*ii-1,3)>0&n_tr3(2*ii,3)>0&minnumerror4(ii)==0
                   n_tr5_1(2*ii-1,1)=sum(n_tr5_1(1:2*(ii-1),3))+1;
                   n_tr5_1(2*ii-1,2)=sum(n_tr5_1(1:2*(ii-1),3))+floor(n_tr3(2*ii-1,3)/2);
                   n_tr5_1(2*ii-1,3)=floor(n_tr3(2*ii-1,3)/2);
                   n_tr5_1(2*ii,1)=sum(n_tr5_1(1:2*ii-1,3))+1;
                   n_tr5_1(2*ii,2)=sum(n_tr5_1(1:2*ii-1,3))+floor(n_tr3(2*ii,3)/2);
                   n_tr5_1(2*ii,3)=floor(n_tr3(2*ii,3)/2);
                   a5_1(n_tr5_1(2*ii-1,1):n_tr5_1(2*ii-1,2),:)=a3(sum(n_tr3(1:2*(ii-1),3))+1:sum(n_tr3(1:2*(ii-1),3))+floor(n_tr3(2*ii-1,3)/2),:);
                   d5_1(n_tr5_1(2*ii-1,1):n_tr5_1(2*ii-1,2),1)=ones(floor(n_tr3(2*ii-1,3)/2),1);
                   a5_1(n_tr5_1(2*ii,1):n_tr5_1(2*ii,2),:)=a3(sum(n_tr3(1:2*ii-1,3))+1:sum(n_tr3(1:2*ii-1,3))+floor(n_tr3(2*ii,3)/2),:);
                   d5_1(n_tr5_1(2*ii,1):n_tr5_1(2*ii,2),1)=-ones(floor(n_tr3(2*ii,3)/2),1);

                   n_tr5_2(2*ii-1,1)=sum(n_tr5_2(1:2*(ii-1),3))+1;
                   n_tr5_2(2*ii-1,2)=sum(n_tr5_2(1:2*(ii-1),3))+n_tr3(2*ii-1,3)-floor(n_tr3(2*ii-1,3)/2);
                   n_tr5_2(2*ii-1,3)=n_tr3(2*ii-1,3)-floor(n_tr3(2*ii-1,3)/2);
                   n_tr5_2(2*ii,1)=sum(n_tr5_2(1:2*ii-1,3))+1;
                   n_tr5_2(2*ii,2)=sum(n_tr5_2(1:2*ii-1,3))+n_tr3(2*ii,3)-floor(n_tr3(2*ii,3)/2);
                   n_tr5_2(2*ii,3)=n_tr3(2*ii,3)-floor(n_tr3(2*ii,3)/2);
                   a5_2(n_tr5_2(2*ii-1,1):n_tr5_2(2*ii-1,2),:)=a3(sum(n_tr3(1:2*(ii-1),3))+floor(n_tr3(2*ii-1,3)/2)+1:n_tr3(2*ii-1,2),:);
                   d5_2(n_tr5_2(2*ii-1,1):n_tr5_2(2*ii-1,2),1)=ones(n_tr3(2*ii-1,3)-floor(n_tr3(2*ii-1,3)/2),1);
                   a5_2(n_tr5_2(2*ii,1):n_tr5_2(2*ii,2),:)=a3(sum(n_tr3(1:2*ii-1,3))+floor(n_tr3(2*ii,3)/2)+1:n_tr3(2*ii,2),:);
                   d5_2(n_tr5_2(2*ii,1):n_tr5_2(2*ii,2),1)=-ones(n_tr3(2*ii,3)-floor(n_tr3(2*ii,3)/2),1);
                end

                if n_tr3(2*ii-1,3)>0&n_tr3(2*ii,3)>0&minnumerror4(ii)>0
                   n_tr4=[];
                   n_tr4(1,1)=1;
                   n_tr4(1,2)=n_tr3(2*ii-1,3);
                   n_tr4(1,3)=n_tr3(2*ii-1,3);
                   n_tr4(1,4)=1;
                   n_tr4(2,1)=n_tr3(2*ii-1,3)+1;
                   n_tr4(2,2)=n_tr3(2*ii-1,3)+n_tr3(2*ii,3);
                   n_tr4(2,3)=n_tr3(2*ii,3);
                   n_tr4(2,4)=2;
                   a00=a3(n_tr3(2*ii-1,1):n_tr3(2*ii,2),:);

                  [a51,d51,n_tr51]=sample_decomposition_7(a00,n_tr4,w5,theta5,dim);
                   if n_tr51(1,3)==0&n_tr51(2,3)==0
                      n_tr5_1(2*ii-1:2*ii,1:3)=zeros(2,3);
                   end
                   if n_tr51(1,3)>0&n_tr51(2,3)==0
                      n_tr5_1(2*ii-1,1:2)=sum(n_tr5_1(1:2*(ii-1),3))+n_tr51(1,1:2);
                      n_tr5_1(2*ii-1,3)=n_tr51(1,3);
                      n_tr5_1(2*ii,1:3)=zeros(1,3);
                      a5_1(n_tr5_1(2*ii-1,1):n_tr5_1(2*ii-1,2),:)=a51(n_tr51(1,1):n_tr51(1,2),:);
                      d5_1(n_tr5_1(2*ii-1,1):n_tr5_1(2*ii-1,2),1)=d51(n_tr51(1,1):n_tr51(1,2),:);
                   end
                   if n_tr51(1,3)==0&n_tr51(2,3)>0
                      n_tr5_1(2*ii-1,1:3)=zeros(1,3);
                      n_tr5_1(2*ii,1:2)=sum(n_tr5_1(1:2*ii-1,3))+n_tr51(2,1:2);
                      n_tr5_1(2*ii,3)=n_tr51(2,3);
                      a5_1(n_tr5_1(2*ii,1):n_tr5_1(2*ii,2),:)=a51(n_tr51(2,1):n_tr51(2,2),:);
                      d5_1(n_tr5_1(2*ii,1):n_tr5_1(2*ii,2),1)=d51(n_tr51(2,1):n_tr51(2,2),:);
                   end
                   if n_tr51(1,3)>0&n_tr51(2,3)>0
                      n_tr5_1(2*ii-1:2*ii,1:2)=sum(n_tr5_1(1:2*(ii-1),3))+n_tr51(1:2,1:2);
                      n_tr5_1(2*ii-1:2*ii,3)=n_tr51(1:2,3);
                      a5_1(n_tr5_1(2*ii-1,1):n_tr5_1(2*ii,2),:)=a51(n_tr51(1,1):n_tr51(2,2),:);
                      d5_1(n_tr5_1(2*ii-1,1):n_tr5_1(2*ii,2),1)=d51(n_tr51(1,1):n_tr51(2,2),:);
                   end

                   if n_tr51(3,3)==0&n_tr51(4,3)==0
                      n_tr5_2(2*ii-1:2*ii,1:3)=zeros(2,3);
                   end
                   if n_tr51(3,3)>0&n_tr51(4,3)==0
                      n_tr5_2(2*ii-1,1)=sum(n_tr5_2(1:2*(ii-1),3))+1;
                      n_tr5_2(2*ii-1,2)=sum(n_tr5_2(1:2*(ii-1),3))+n_tr51(3,3);
                      n_tr5_2(2*ii-1,3)=n_tr51(3,3);
                      n_tr5_2(2*ii,1:3)=zeros(1,3);
                      a5_2(n_tr5_2(2*ii-1,1):n_tr5_2(2*ii-1,2),:)=a51(n_tr51(3,1):n_tr51(3,2),:);
                      d5_2(n_tr5_2(2*ii-1,1):n_tr5_2(2*ii-1,2),1)=d51(n_tr51(3,1):n_tr51(3,2),:);
                   end
                   if n_tr51(3,3)==0&n_tr51(4,3)>0
                      n_tr5_2(2*ii-1,1:3)=zeros(1,3);
                      n_tr5_2(2*ii,1)=sum(n_tr5_2(1:2*ii-1,3))+1;
                      n_tr5_2(2*ii,2)=sum(n_tr5_2(1:2*ii-1,3))+n_tr51(4,3);
                      n_tr5_2(2*ii,3)=n_tr51(4,3);
                      a5_2(n_tr5_2(2*ii,1):n_tr5_2(2*ii,2),:)=a51(n_tr51(4,1):n_tr51(4,2),:);
                      d5_2(n_tr5_2(2*ii,1):n_tr5_2(2*ii,2),1)=d51(n_tr51(4,1):n_tr51(4,2),:);
                   end
                   if n_tr51(3,3)>0&n_tr51(4,3)>0
                      n_tr5_2(2*ii-1,1)=sum(n_tr5_2(1:2*(ii-1),3))+1;
                      n_tr5_2(2*ii-1,2)=sum(n_tr5_2(1:2*(ii-1),3))+n_tr51(3,3);
                      n_tr5_2(2*ii-1:2*ii,3)=n_tr51(3:4,3);
                      n_tr5_2(2*ii,1)=sum(n_tr5_2(1:2*ii-1,3))+1;
                      n_tr5_2(2*ii,2)=sum(n_tr5_2(1:2*ii-1,3))+n_tr51(4,3);
                      a5_2(n_tr5_2(2*ii-1,1):n_tr5_2(2*ii,2),:)=a51(n_tr51(3,1):n_tr51(4,2),:);
                      d5_2(n_tr5_2(2*ii-1,1):n_tr5_2(2*ii,2),1)=d51(n_tr51(3,1):n_tr51(4,2),:);
                   end
                end
            end
            n_tr5=n_tr5_1;
            n_tr5(9:16,1:2)=ones(8,2)*sum(n_tr5_1(1:8,3))+n_tr5_2(:,1:2);
            n_tr5(9:16,3)=n_tr5_2(:,3);
            for ii=1:16,
                n_tr5(ii,4)=ii;
            end
            for ii=1:16,
                if n_tr5(ii,3)==0
                   n_tr5(ii,1:2)=zeros(1,2);
                end
            end       

            n_tr5
            fprintf(fid,'\n n_tr5=\n',n_tr5');
            fprintf(fid,'%g  %g  %g  %g\n',n_tr5');
            fprintf(fid,'\n');

            a5=a5_1;
            a5(sum(n_tr5(1:8,3))+1:sum(n_tr5(:,3)),:)=a5_2;
            d5=d5_1;
            d5(sum(n_tr5(1:8,3))+1:sum(n_tr5(:,3)),:)=d5_2;
            minnumerror5_1=zeros(1,8);
            minnumerror5_2=zeros(8,2);
            for i=1:2,
                if n_tr5(8*i-7,3)==0&n_tr5(8*i-6,3)>0
                   w5_1(4*i-3,:)=-w2(k0,:);
                   theta5_1(4*i-3)=-theta2(k0);
                   n0000(4*i-3)=0;
                end
                if n_tr5(8*i-7,3)>0&n_tr5(8*i-6,3)==0
                   w5_1(4*i-3,:)=w2(k0,:);
                   theta5_1(4*i-3)=theta2(k0);
                   n0000(4*i-3)=0;
                end
                if n_tr5(8*i-5,3)==0&n_tr5(8*i-4,3)>0
                   w5_1(4*i-2,:)=w2(k0,:);
                   theta5_1(4*i-2)=theta2(k0);
                   n0000(4*i-2)=0;
                end
                if n_tr5(8*i-5,3)>0&n_tr5(8*i-4,3)==0
                   w5_1(4*i-2,:)=-w2(k0,:);
                   theta5_1(4*i-2)=-theta2(k0);
                   n0000(4*i-2)=0;
                end
                if n_tr5(8*i-3,3)==0&n_tr5(8*i-2,3)>0
                   w5_1(4*i-1,:)=-w2(k0,:);
                   theta5_1(4*i-1)=mu01*w2(k0,:)';
                   n0000(4*i-1)=0;
                end
                if n_tr5(8*i-3,3)>0&n_tr5(8*i-2,3)==0
                   w5_1(4*i-1,:)=w2(k0,:);
                   theta5_1(4*i-1)=-mu01*w2(k0,:)';
                   n0000(4*i-1)=0;
                end

                if n_tr5(8*i-1,3)==0&n_tr5(8*i,3)>0
                   w5_1(4*i,:)=w2(k0,:);
                   theta5_1(4*i)=-mu02*w2(k0,:)';
                   n0000(4*i)=0;
                end
                if n_tr5(8*i-1,3)>0&n_tr5(8*i,3)==0
                   w5_1(4*i,:)=-w2(k0,:);
                   theta5_1(4*i)=mu02*w2(k0,:)';
                   n0000(4*i)=0;
                end
            end

            for ii=1:8,
                p2=n_tr5(2*ii-1,3);
                q2=n_tr5(2*ii,3);
                a22=[];
                d22=[];
                if p2>0&q2>0&((ii<=4&minnumerror4(ii)>0)|(ii>4&minnumerror4(ii-4)>0))
                   a22=a5(n_tr5(2*ii-1,1):n_tr5(2*ii,2),:);
                   d22=d5(n_tr5(2*ii-1,1):n_tr5(2*ii,2),:);
                  [w3,theta3,minnumerror3,n]=calculating_weigh_theta_6(a22,dim,d22,p0,q0,p2,q2,fid);
%                   minnumerror3_0=floor(10*minnumerror3(1));
                   fprintf('num(%d,%d,%d)=%d\n\n',r,t,ii,n);
                   fprintf(fid,'num(%d,%d,%d)=%d\n\n',r,t,ii,n);
                   w5_1(ii,:)=w3';
                   theta5_1(ii)=theta3;
                   minnumerror5_1(ii)=sum(minnumerror3(2:3));
                   minnumerror5_2(ii,:)=minnumerror3(2:3);
                   n0000(ii)=n;
                end
            end
            n000(k0)=n000(k0)+sum(n0000);

            minnumerror5=minnumerror5_1;
            fprintf('minnumerror5=%g  %g  %g  %g  %g  %g  %g  %g, sum=%g\n',minnumerror5,sum(minnumerror5));
            fprintf(fid,'minnumerror5=%g  %g  %g  %g  %g  %g  %g  %g, sum=%g\n',minnumerror5,sum(minnumerror5));

            s1=sum(sign0)-6;
            w6(s1,:)=w2(k0,:);
            theta6(s1)=theta2(k0);

                sign0(k0)=12;
                minnumerror6_1=[0 0 0 0 0 0 0 0];

                w6_1=w5_1;
                theta6_1=theta5_1;
                minnumerror6_1=minnumerror5_1;
                 
                if minnumerror5_2(1,1)+minnumerror5_2(1,2)>=n_tr5(2,3)
                   w6_1(1,:)=w44(1,:);
                   theta6_1(1)=theta44(1);
                   minnumerror6_1(1)=n_tr5(2,3);
                end
                if minnumerror5_2(5,1)+minnumerror5_2(5,2)>=n_tr5(10,3)
                   w6_1(5,:)=w44(1,:);
                   theta6_1(5)=theta44(1);
                   minnumerror6_1(5)=n_tr5(10,3);
                end
                if minnumerror4(1)==0
                    w6_1(1,:)=w44(1,:);
                    theta6_1(1)=theta44(1);
                    w6_1(5,:)=w44(1,:);
                    theta6_1(5)=theta44(1);
                    minnumerror6_1(1)=0;
                    minnumerror6_1(5)=0;
                end

                if minnumerror5_2(2,1)+minnumerror5_2(2,2)>=n_tr5(3,3)
                    w6_1(2,:)=w44(2,:);
                    theta6_1(2)=theta44(2);
                    minnumerror6_1(2)=n_tr5(3,3);
                end
                if minnumerror5_2(6,1)+minnumerror5_2(6,2)>=n_tr5(11,3)
                    w6_1(6,:)=w44(2,:);
                    theta6_1(6)=theta44(2);
                    minnumerror6_1(6)=n_tr5(11,3);
                end
                if minnumerror4(2)==0
                    w6_1(2,:)=w44(2,:);
                    theta6_1(2)=theta44(2);
                    w6_1(6,:)=w44(2,:);
                    theta6_1(6)=theta44(2);
                    minnumerror6_1(2)=0;
                    minnumerror6_1(6)=0;
                end

                if minnumerror5_2(3,1)+minnumerror5_2(3,2)>=n_tr5(6,3)
                    w6_1(3,:)=w44(3,:);
                    theta6_1(3)=theta44(3);
                    minnumerror6_1(3)=n_tr5(6,3);
                end
                if minnumerror5_2(7,1)+minnumerror5_2(7,2)>=n_tr5(14,3)
                    w6_1(7,:)=w44(3,:);
                    theta6_1(7)=theta44(3);
                    minnumerror6_1(7)=n_tr5(14,3);
                end
                if minnumerror4(3)==0
                    w6_1(3,:)=w44(3,:);
                    theta6_1(3)=theta44(3);
                    w6_1(7,:)=w44(3,:);
                    theta6_1(7)=theta44(3);
                    minnumerror6_1(3)=0;
                    minnumerror6_1(7)=0;
                end

                if minnumerror5_2(4,1)+minnumerror5_2(4,2)>=n_tr5(7,3)
                    w6_1(4,:)=w44(4,:);
                    theta6_1(4)=theta44(4);
                    minnumerror6_1(4)=n_tr5(7,3);
                end
                if minnumerror5_2(8,1)+minnumerror5_2(8,2)>=n_tr5(15,3)
                    w6_1(8,:)=w44(4,:);
                    theta6_1(8)=theta44(4);
                    minnumerror6_1(8)=n_tr5(15,3);
                end
                if minnumerror4(4)==0
                    w6_1(4,:)=w44(4,:);
                    theta6_1(4)=theta44(4);
                    w6_1(8,:)=w44(4,:);
                    theta6_1(8)=theta44(4);
                    minnumerror6_1(4)=0;
                    minnumerror6_1(8)=0;
                end

                w6(s1+4:s1+11,:)=w6_1;
                theta6(s1+4:s1+11)=theta6_1;
                minnumerror6=minnumerror6_1;
                w66=w6_1;
                theta66=theta6_1;
            w6(s1+1,:)=w2(k0,:);
            theta6(s1+1)=-mu01*w2(k0,:)';
            w6(s1+2,:)=w2(k0,:);
            theta6(s1+2)=-mu02*w2(k0,:)';
            w6(s1+3,:)=w5;
            theta6(s1+3)=theta5;
            fprintf('minnumerror6=%g  %g  %g  %g  %g  %g  %g  %g, sum=%g\n',minnumerror6,sum(minnumerror6));
            fprintf(fid,'minnumerror6=%g  %g  %g  %g  %g  %g  %g  %g, sum=%g\n',minnumerror6,sum(minnumerror6));

            kk=sum(sign0(1:k0));
            if sign0(k0)==12
               numtr0=[0 0];
               re1=a0(:,2:dim+1)*(w6(kk-11:kk,:))'+ones(p0+q0,1)*theta6(kk-11:kk);
               numtr1=zeros(8,2);
               for i=1:p0,
                   if (re1(i,4)>0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,5)<0)
                       numtr0(1)=numtr0(1)+1;
                       numtr1(1,1)=numtr1(1,1)+1;
                   end
                   if (re1(i,4)>0)&(re1(i,1)<0)&(re1(i,3)>0)&(re1(i,6)<0)
                       numtr0(1)=numtr0(1)+1;
                       numtr1(2,1)=numtr1(2,1)+1;
                   end
                   if (re1(i,4)>0)&(re1(i,1)>0)&(re1(i,2)>0)&(re1(i,7)<0)
                       numtr0(1)=numtr0(1)+1;
                       numtr1(3,1)=numtr1(3,1)+1;
                   end
                   if (re1(i,4)>0)&(re1(i,1)<0)&(re1(i,3)<0)&(re1(i,8)<0)
                       numtr0(1)=numtr0(1)+1;
                       numtr1(4,1)=numtr1(4,1)+1;
                   end
                   if (re1(i,4)<0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,9)<0)
                       numtr0(1)=numtr0(1)+1;
                       numtr1(5,1)=numtr1(5,1)+1;
                   end
                   if (re1(i,4)<0)&(re1(i,1)<0)&(re1(i,3)>0)&(re1(i,10)<0)
                       numtr0(1)=numtr0(1)+1;
                       numtr1(6,1)=numtr1(6,1)+1;
                   end
                   if (re1(i,4)<0)&(re1(i,1)>0)&(re1(i,2)>0)&(re1(i,11)<0)
                       numtr0(1)=numtr0(1)+1;
                       numtr1(7,1)=numtr1(7,1)+1;
                   end
                   if (re1(i,4)<0)&(re1(i,1)<0)&(re1(i,3)<0)&(re1(i,12)<0)
                       numtr0(1)=numtr0(1)+1;
                       numtr1(8,1)=numtr1(8,1)+1;
                   end
               end

               for i=p0+1:p0+q0,
                   if (re1(i,4)>0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,5)>0)
                       numtr0(2)=numtr0(2)+1;
                       numtr1(1,2)=numtr1(1,2)+1;
                   end
                   if (re1(i,4)>0)&(re1(i,1)<0)&(re1(i,3)>0)&(re1(i,6)>0)
                       numtr0(2)=numtr0(2)+1;
                       numtr1(2,2)=numtr1(2,2)+1;
                   end
                   if (re1(i,4)>0)&(re1(i,1)>0)&(re1(i,2)>0)&(re1(i,7)>0)
                       numtr0(2)=numtr0(2)+1;
                       numtr1(3,2)=numtr1(3,2)+1;
                   end
                   if (re1(i,4)>0)&(re1(i,1)<0)&(re1(i,3)<0)&(re1(i,8)>0)
                       numtr0(2)=numtr0(2)+1;
                       numtr1(4,2)=numtr1(4,2)+1;
                   end
                   if (re1(i,4)<0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,9)>0)
                       numtr0(2)=numtr0(2)+1;
                       numtr1(5,2)=numtr1(5,2)+1;
                   end
                   if (re1(i,4)<0)&(re1(i,1)<0)&(re1(i,3)>0)&(re1(i,10)>0)
                       numtr0(2)=numtr0(2)+1;
                       numtr1(6,2)=numtr1(6,2)+1;
                   end
                   if (re1(i,4)<0)&(re1(i,1)>0)&(re1(i,2)>0)&(re1(i,11)>0)
                       numtr0(2)=numtr0(2)+1;
                       numtr1(7,2)=numtr1(7,2)+1;
                   end
                   if (re1(i,4)<0)&(re1(i,1)<0)&(re1(i,3)<0)&(re1(i,12)>0)
                       numtr0(2)=numtr0(2)+1;
                       numtr1(8,2)=numtr1(8,2)+1;
                   end
               end
               numtr0
               numtr1'
               sum(sum(numtr1))
               minnumerror6=sum(numtr1');
               fprintf('minnumerror6=%g  %g  %g  %g  %g  %g  %g  %g, sum=%g\n',minnumerror6,sum(minnumerror6));
               fprintf(fid,'minnumerror6=%g  %g  %g  %g  %g  %g  %g  %g, sum=%g\n',minnumerror6,sum(minnumerror6));
               minnumerror2(k0)=sum(numtr0);
               fprintf(fid,'%g  %g\n',numtr0');
               fprintf(fid,'%g  %g\n',numtr1');
               fprintf(fid,'%g\n',sum(sum(numtr1)));
            end
         end
         numtr_0(3)=numtr0(1);
         tacc_index(3)=sum(numtr0);
%%////////////////////////////////////////////////////////////////////////////
        pe=numtr0(1)/p0;
        qe=numtr0(2)/q0;
        if ir<=19
           aver_index(3)=(numtr0(1)+numtr0(2))/(p0+q0);
        elseif ir<=49
           aver_index(3)=(pe+qe)/2;
        elseif ir<=99.00
           if numtr0(1)>0&numtr0(2)>0
              aver_index(3)=(1-sqrt((1-pe)*(1-qe)));
           elseif numtr0(1)>0&numtr0(2)==0
              aver_index(3)=(1-sqrt(0.0001*(1-pe)));
           elseif numtr0(1)==0&numtr0(2)>0
              aver_index(3)=(1-sqrt(0.0001*(1-qe)));
           else
              aver_index(3)=0;
           end
        else
           if numtr0(1)>0&numtr0(2)>0
              aver_index(3)=(1/(1+alpha))*((pe+beta*qe-2*beta*pe*qe)/(2-pe-beta*qe)+alpha*pe);
           elseif numtr0(1)>0&numtr0(2)==0
              aver_index(3)=(1/(1+alpha))*((pe+beta*0.0001-2*beta*pe*0.0001)/(2-pe-beta*0.0001)+alpha*pe);
           elseif numtr0(1)==0&numtr0(2)>0
              aver_index(3)=(1/(1+alpha))*((0.0001+beta*qe-2*beta*0.0001*qe)/(2-0.0001-beta*qe)+alpha*pe);
           else
              aver_index(3)=0;
           end
        end
%%/////////////////////////////////////////////////////////////////////////////////////////////
         t1=etime(clock,t0)
         fprintf(fid,'runtime(%d,%d)=%g\n',r,t,t1);
%///////////////////////////////////////////////////////////////////////////////////
         if (aver_index(3)>0.05)&(((ir>99)&(numtr_0(3)<=numtr_0(2))&(numtr_0(2)<=numtr_0(1)))|(ir<=99))

            project3=a5(:,2:dim+1)*w5';
            xigma3=std(project3);

            a7=[];
            d7=[];
            n_tr7=[];
            a7_1=[];
            a7_2=[];
            d7_1=[];
            d7_2=[];
            n_tr7_1=[];
            n_tr7_2=[];

            n_tr7_1=zeros(8,3);
            n_tr7_2=zeros(8,3);
            if n_tr5(1,3)==0&n_tr5(2,3)>0
               n_tr7_1(1,1:3)=zeros(1,3);
               n_tr7_1(2,1)=1;
               n_tr7_1(2,2)=floor(n_tr5(2,3)/2);
               n_tr7_1(2,3)=n_tr7_1(2,2);
               a7_1=a5(n_tr5(2,1):floor(n_tr5(2,3)/2),:);
               d7_1(:,1)=-ones(floor(n_tr5(2,3)/2),1);

               n_tr7_2(1,1:3)=zeros(1,3);
               n_tr7_2(2,1)=1;
               n_tr7_2(2,2)=n_tr5(2,3)-floor(n_tr5(2,3)/2);
               n_tr7_2(2,3)=n_tr7_2(2,2);
               a7_2=a5(floor(n_tr5(2,3)/2)+1:n_tr5(2,2),:);
               d7_2(:,1)=-ones(n_tr5(2,3)-floor(n_tr5(2,3)/2),1);
            end
            if n_tr5(1,3)>0&n_tr5(2,3)==0
               if n_tr5(1,3)>1
                  n_tr7_1(1,1)=1;
                  n_tr7_1(1,2)=floor(n_tr5(1,3)/2);
                  n_tr7_1(1,3)=n_tr7_1(1,2);
                  a7_1=a5(n_tr5(1,1):floor(n_tr5(1,3)/2),:);
                  d7_1(:,1)=ones(floor(n_tr5(1,3)/2),1);
               else
                  n_tr7_1(1,1:3)=zeros(1,3);
               end
               n_tr7_1(2,1:3)=zeros(1,3);
               n_tr7_2(1,1)=1;
               n_tr7_2(1,2)=n_tr5(1,3)-floor(n_tr5(1,3)/2);
               n_tr7_2(1,3)=n_tr7_2(1,2);
               n_tr7_2(2,1:3)=zeros(1,3);
               a7_2=a5(floor(n_tr5(1,3)/2)+1:n_tr5(1,3),:);
               d7_2(:,1)=ones(n_tr5(1,3)-floor(n_tr5(1,3)/2),1);
            end

            if (n_tr5(1,3)>0&n_tr5(2,3)>0)&minnumerror6(1)==0
               if n_tr5(1,3)>1
                  n_tr7_1(1,1)=1;
                  n_tr7_1(1,2)=floor(n_tr5(1,3)/2);
                  n_tr7_1(1,3)=n_tr7_1(1,2);
                  a7_1=a5(n_tr5(1,1):floor(n_tr5(1,3)/2),:);
                  d7_1=ones(floor(n_tr5(1,3)/2),1);
               else
                  n_tr7_1(1,1:3)=zeros(1,3);
               end
               n_tr7_1(2,1)=n_tr7_1(1,2)+1;
               n_tr7_1(2,2)=n_tr7_1(1,2)+floor(n_tr5(2,3)/2);
               n_tr7_1(2,3)=floor(n_tr5(2,3)/2);

               a7_1(n_tr7_1(2,1):n_tr7_1(2,2),:)=a5(n_tr5(2,1):n_tr5(2,1)+floor(n_tr5(2,3)/2)-1,:);
               d7_1(n_tr7_1(2,1):n_tr7_1(2,2),1)=-ones(floor(n_tr5(2,3)/2),1);

               n_tr7_2(1,1)=1;
               n_tr7_2(1,2)=n_tr5(1,3)-floor(n_tr5(1,3)/2);
               n_tr7_2(1,3)=n_tr7_2(1,2);
               n_tr7_2(2,1)=n_tr7_2(1,2)+1;
               n_tr7_2(2,2)=n_tr7_2(1,2)+n_tr5(2,3)-floor(n_tr5(2,3)/2);
               n_tr7_2(2,3)=n_tr5(2,3)-floor(n_tr5(2,3)/2);

               a7_2=a5(floor(n_tr5(1,3)/2)+1:n_tr5(1,3),:);
               a7_2(n_tr7_2(2,1):n_tr7_2(2,2),:)=a5(n_tr5(1,2)+floor(n_tr5(2,3)/2)+1:n_tr5(2,2),:);
               d7_2(:,1)=ones(n_tr5(1,3)-floor(n_tr5(1,3)/2),1);
               d7_2(n_tr7_2(2,1):n_tr7_2(2,2),1)=-ones(n_tr5(2,3)-floor(n_tr5(2,3)/2),1);
            end

            if (n_tr5(1,3)>0&n_tr5(2,3)>0)&minnumerror6(1)>0
               n_tr6=[];
               n_tr6(1,1)=1;
               n_tr6(1,2)=n_tr5(1,3);
               n_tr6(1,3)=n_tr5(1,3);
               n_tr6(1,4)=1;
               n_tr6(2,1)=n_tr5(1,3)+1;
               n_tr6(2,2)=n_tr5(1,3)+n_tr5(2,3);
               n_tr6(2,3)=n_tr5(2,3);
               n_tr6(2,4)=2;
               a00=[];
               a00=a5(n_tr5(1,1):n_tr5(2,2),:);

              [a71,d71,n_tr71]=sample_decomposition_8(a00,n_tr6,w5,theta5,xigma3,dim);
               n_tr7_1(1:2,1:3)=n_tr71(1:2,1:3);

               if n_tr71(1,3)==0&n_tr71(2,3)>0
                  a7_1=a71(n_tr71(2,1):n_tr71(2,2),:);
                  d7_1=d71(n_tr71(2,1):n_tr71(2,2),:);
               end
               if n_tr71(1,3)>0&n_tr71(2,3)==0
                  a7_1=a71(n_tr71(1,1):n_tr71(1,2),:);
                  d7_1=d71(n_tr71(1,1):n_tr71(1,2),:);
               end
               if n_tr71(1,3)>0&n_tr71(2,3)>0
                  a7_1=a71(n_tr71(1,1):n_tr71(2,2),:);
                  d7_1=d71(n_tr71(1,1):n_tr71(2,2),:);
               end

               if n_tr71(3,3)==0&n_tr71(4,3)>0
                  n_tr7_2(1,1:3)=zeros(1,3);
                  n_tr7_2(2,1)=1;
                  n_tr7_2(2,2)=n_tr71(4,3);
                  n_tr7_2(2,3)=n_tr71(4,3);
                  a7_2=a71(n_tr71(4,1):n_tr71(4,2),:);
                  d7_2=d71(n_tr71(4,1):n_tr71(4,2),:);
               end
               if n_tr71(3,3)>0&n_tr71(4,3)==0
                  n_tr7_2(1,1)=1;
                  n_tr7_2(1,2)=n_tr71(3,3);
                  n_tr7_2(1,3)=n_tr71(3,3);
                  n_tr7_2(2,1:3)=zeros(1,3);
                  a7_2=a71(n_tr71(3,1):n_tr71(3,2),:);
                  d7_2=d71(n_tr71(3,1):n_tr71(3,2),:);
               end
               if n_tr71(3,3)>0&n_tr71(4,3)>0
                  n_tr7_2(1,1)=1;
                  n_tr7_2(1,2)=n_tr71(3,3);
                  n_tr7_2(2,1)=n_tr71(3,3)+1;
                  n_tr7_2(2,2)=n_tr71(3,3)+n_tr71(4,3);
                  n_tr7_2(1:2,3)=n_tr71(3:4,3);
                  a7_2=a71(n_tr71(3,1):n_tr71(4,2),:);
                  d7_2=d71(n_tr71(3,1):n_tr71(4,2),:);
               end
            end

            for ii=2:4,
                n_tr7_1(2*ii-1:2*ii,1:3)=zeros(2,3);
                n_tr7_2(2*ii-1:2*ii,1:3)=zeros(2,3);
                if n_tr5(2*ii-1,3)==0&n_tr5(2*ii,3)==0
                   n_tr7_1(2*ii-1:2*ii,1:3)=zeros(2,3);
                   n_tr7_2(2*ii-1:2*ii,1:3)=zeros(2,3);
                end
                if n_tr5(2*ii-1,3)==0&n_tr5(2*ii,3)>0
                   n_tr7_1(2*ii-1,1:3)=zeros(1,3);
                   n_tr7_1(2*ii,1)=sum(n_tr7_1(1:2*(ii-1),3))+1;
                   n_tr7_1(2*ii,2)=sum(n_tr7_1(1:2*(ii-1),3))+floor(n_tr5(2*ii,3)/2);
                   n_tr7_1(2*ii,3)=floor(n_tr5(2*ii,3)/2);
                   a7_1(n_tr7_1(2*ii,1):n_tr7_1(2*ii,2),:)=a5(sum(n_tr5(1:2*(ii-1),3))+1:sum(n_tr5(1:2*(ii-1),3))+floor(n_tr5(2*ii,3)/2),:);
                   d7_1(n_tr7_1(2*ii,1):n_tr7_1(2*ii,2),1)=-ones(floor(n_tr5(2*ii,3)/2),1);

                   n_tr7_2(2*ii-1,1:3)=zeros(1,3);
                   n_tr7_2(2*ii,1)=sum(n_tr7_2(1:2*(ii-1),3))+1;
                   n_tr7_2(2*ii,2)=sum(n_tr7_2(1:2*(ii-1),3))+n_tr5(2*ii,3)-floor(n_tr5(2*ii,3)/2);
                   n_tr7_2(2*ii,3)=n_tr5(2*ii,3)-floor(n_tr5(2*ii,3)/2);
                   a7_2(n_tr7_2(2*ii,1):n_tr7_2(2*ii,2),:)=a5(sum(n_tr5(1:2*(ii-1),3))+floor(n_tr5(2*ii,3)/2)+1:sum(n_tr5(1:2*(ii-1),3))+n_tr5(2*ii,3),:);
                   d7_2(n_tr7_2(2*ii,1):n_tr7_2(2*ii,2),1)=-ones(n_tr5(2*ii,3)-floor(n_tr5(2*ii,3)/2),1);
                end
                if n_tr5(2*ii-1,3)>0&n_tr5(2*ii,3)==0
                   n_tr7_1(2*ii-1,1)=sum(n_tr7_1(1:2*(ii-1),3))+1;
                   n_tr7_1(2*ii-1,2)=sum(n_tr7_1(1:2*(ii-1),3))+floor(n_tr5(2*ii-1,3)/2);
                   n_tr7_1(2*ii-1,3)=floor(n_tr5(2*ii-1,3)/2);
                   a7_1(n_tr7_1(2*ii-1,1):n_tr7_1(2*ii-1,2),:)=a5(sum(n_tr5(1:2*(ii-1),3))+1:sum(n_tr5(1:2*(ii-1),3))+floor(n_tr5(2*ii-1,3)/2),:);
                   d7_1(n_tr7_1(2*ii-1,1):n_tr7_1(2*ii-1,2),1)=ones(floor(n_tr5(2*ii-1,3)/2),1);
                   n_tr7_1(2*ii,1:3)=zeros(1,3);

                   n_tr7_2(2*ii-1,1)=sum(n_tr7_2(1:2*(ii-1),3))+1;
                   n_tr7_2(2*ii-1,2)=sum(n_tr7_2(1:2*(ii-1),3))+n_tr5(2*ii-1,3)-floor(n_tr5(2*ii-1,3)/2);
                   n_tr7_2(2*ii-1,3)=n_tr5(2*ii-1,3)-floor(n_tr5(2*ii-1,3)/2);
                   a7_2(n_tr7_2(2*ii-1,1):n_tr7_2(2*ii-1,2),:)=a5(sum(n_tr5(1:2*(ii-1),3))+floor(n_tr5(2*ii-1,3)/2)+1:sum(n_tr5(1:2*(ii-1),3))+n_tr5(2*ii-1,3),:);
                   d7_2(n_tr7_2(2*ii-1,1):n_tr7_2(2*ii-1,2),1)=ones(n_tr5(2*ii-1,3)-floor(n_tr5(2*ii-1,3)/2),1);
                   n_tr7_2(2*ii,1:3)=zeros(1,3);
                end

                if n_tr5(2*ii-1,3)>0&n_tr5(2*ii,3)>0&minnumerror6(ii)==0
                   n_tr7_1(2*ii-1,1)=sum(n_tr7_1(1:2*(ii-1),3))+1;
                   n_tr7_1(2*ii-1,2)=sum(n_tr7_1(1:2*(ii-1),3))+floor(n_tr5(2*ii-1,3)/2);
                   n_tr7_1(2*ii-1,3)=floor(n_tr5(2*ii-1,3)/2);
                   n_tr7_1(2*ii,1)=sum(n_tr7_1(1:2*ii-1,3))+1;
                   n_tr7_1(2*ii,2)=sum(n_tr7_1(1:2*ii-1,3))+floor(n_tr5(2*ii,3)/2);
                   n_tr7_1(2*ii,3)=floor(n_tr5(2*ii,3)/2);
                   a7_1(n_tr7_1(2*ii-1,1):n_tr7_1(2*ii-1,2),:)=a5(sum(n_tr5(1:2*(ii-1),3))+1:sum(n_tr5(1:2*(ii-1),3))+floor(n_tr5(2*ii-1,3)/2),:);
                   d7_1(n_tr7_1(2*ii-1,1):n_tr7_1(2*ii-1,2),1)=ones(floor(n_tr5(2*ii-1,3)/2),1);
                   a7_1(n_tr7_1(2*ii,1):n_tr7_1(2*ii,2),:)=a5(sum(n_tr5(1:2*ii-1,3))+1:sum(n_tr5(1:2*ii-1,3))+floor(n_tr5(2*ii,3)/2),:);
                   d7_1(n_tr7_1(2*ii,1):n_tr7_1(2*ii,2),1)=-ones(floor(n_tr5(2*ii,3)/2),1);

                   n_tr7_2(2*ii-1,1)=sum(n_tr7_2(1:2*(ii-1),3))+1;
                   n_tr7_2(2*ii-1,2)=sum(n_tr7_2(1:2*(ii-1),3))+n_tr5(2*ii-1,3)-floor(n_tr5(2*ii-1,3)/2);
                   n_tr7_2(2*ii-1,3)=n_tr5(2*ii-1,3)-floor(n_tr5(2*ii-1,3)/2);
                   n_tr7_2(2*ii,1)=sum(n_tr7_2(1:2*ii-1,3))+1;
                   n_tr7_2(2*ii,2)=sum(n_tr7_2(1:2*ii-1,3))+n_tr5(2*ii,3)-floor(n_tr5(2*ii,3)/2);
                   n_tr7_2(2*ii,3)=n_tr5(2*ii,3)-floor(n_tr5(2*ii,3)/2);
                   a7_2(n_tr7_2(2*ii-1,1):n_tr7_2(2*ii-1,2),:)=a5(sum(n_tr5(1:2*(ii-1),3))+floor(n_tr5(2*ii-1,3)/2)+1:n_tr5(2*ii-1,2),:);
                   d7_2(n_tr7_2(2*ii-1,1):n_tr7_2(2*ii-1,2),1)=ones(n_tr5(2*ii-1,3)-floor(n_tr5(2*ii-1,3)/2),1);
                   a7_2(n_tr7_2(2*ii,1):n_tr7_2(2*ii,2),:)=a5(sum(n_tr5(1:2*ii-1,3))+floor(n_tr5(2*ii,3)/2)+1:n_tr5(2*ii,2),:);
                   d7_2(n_tr7_2(2*ii,1):n_tr7_2(2*ii,2),1)=-ones(n_tr5(2*ii,3)-floor(n_tr5(2*ii,3)/2),1);
                end

                if n_tr5(2*ii-1,3)>0&n_tr5(2*ii,3)>0&minnumerror6(ii)>0
                   n_tr6=[];
                   n_tr6(1,1)=1;
                   n_tr6(1,2)=n_tr5(2*ii-1,3);
                   n_tr6(1,3)=n_tr5(2*ii-1,3);
                   n_tr6(1,4)=1;
                   n_tr6(2,1)=n_tr5(2*ii-1,3)+1;
                   n_tr6(2,2)=n_tr5(2*ii-1,3)+n_tr5(2*ii,3);
                   n_tr6(2,3)=n_tr5(2*ii,3);
                   n_tr6(2,4)=2;
                   a00=[];
                   a00=a5(n_tr5(2*ii-1,1):n_tr5(2*ii,2),:);
                  [a71,d71,n_tr71]=sample_decomposition_8(a00,n_tr6,w5,theta5,xigma3,dim);
                   n_tr7_1(2*ii-1:2*ii,1:2)=sum(n_tr7_1(1:2*(ii-1),3))+n_tr71(1:2,1:2);
                   n_tr7_1(2*ii-1:2*ii,3)=n_tr71(1:2,3);
                   if n_tr71(1,3)==0&n_tr71(2,3)>0
                      a7_1(n_tr7_1(2*ii,1):n_tr7_1(2*ii,2),:)=a71(n_tr71(2,1):n_tr71(2,2),:);
                      d7_1(n_tr7_1(2*ii,1):n_tr7_1(2*ii,2),1)=d71(n_tr71(2,1):n_tr71(2,2),:);
                   end
                   if n_tr71(1,3)>0&n_tr71(2,3)==0
                      a7_1(n_tr7_1(2*ii-1,1):n_tr7_1(2*ii-1,2),:)=a71(n_tr71(1,1):n_tr71(1,2),:);
                      d7_1(n_tr7_1(2*ii-1,1):n_tr7_1(2*ii-1,2),1)=d71(n_tr71(1,1):n_tr71(1,2),:);
                   end
                   if n_tr71(1,3)>0&n_tr71(2,3)>0
                      a7_1(n_tr7_1(2*ii-1,1):n_tr7_1(2*ii,2),:)=a71(n_tr71(1,1):n_tr71(2,2),:);
                      d7_1(n_tr7_1(2*ii-1,1):n_tr7_1(2*ii,2),1)=d71(n_tr71(1,1):n_tr71(2,2),:);
                   end

                   if n_tr71(3,3)==0&n_tr71(4,3)>0
                      n_tr7_2(2*ii,3)=n_tr71(4,3);
                      n_tr7_2(2*ii,1)=sum(n_tr7_2(1:2*ii-1,3))+1;
                      n_tr7_2(2*ii,2)=sum(n_tr7_2(1:2*ii-1,3))+n_tr71(4,3);
                      a7_2(n_tr7_2(2*ii,1):n_tr7_2(2*ii,2),:)=a71(n_tr71(4,1):n_tr71(4,2),:);
                      d7_2(n_tr7_2(2*ii,1):n_tr7_2(2*ii,2),1)=d71(n_tr71(4,1):n_tr71(4,2),:);
                   end
                   if n_tr71(3,3)>0&n_tr71(4,3)==0
                      n_tr7_2(2*ii-1,1)=sum(n_tr7_2(1:2*(ii-1),3))+1;
                      n_tr7_2(2*ii-1,2)=sum(n_tr7_2(1:2*(ii-1),3))+n_tr71(3,3);
                      n_tr7_2(2*ii-1,3)=n_tr71(3,3);
                      a7_2(n_tr7_2(2*ii-1,1):n_tr7_2(2*ii-1,2),:)=a71(n_tr71(3,1):n_tr71(3,2),:);
                      d7_2(n_tr7_2(2*ii-1,1):n_tr7_2(2*ii-1,2),1)=d71(n_tr71(3,1):n_tr71(3,2),:);
                   end

                   if n_tr71(3,3)>0&n_tr71(4,3)>0
                      n_tr7_2(2*ii-1:2*ii,3)=n_tr71(3:4,3);
                      n_tr7_2(2*ii-1,1)=sum(n_tr7_2(1:2*(ii-1),3))+1;
                      n_tr7_2(2*ii-1,2)=sum(n_tr7_2(1:2*(ii-1),3))+n_tr71(3,3);
                      n_tr7_2(2*ii,1)=sum(n_tr7_2(1:2*ii-1,3))+1;
                      n_tr7_2(2*ii,2)=sum(n_tr7_2(1:2*ii-1,3))+n_tr71(4,3);
                      a7_2(n_tr7_2(2*ii-1,1):n_tr7_2(2*ii,2),:)=a71(n_tr71(3,1):n_tr71(4,2),:);
                      d7_2(n_tr7_2(2*ii-1,1):n_tr7_2(2*ii,2),1)=d71(n_tr71(3,1):n_tr71(4,2),:);
                   end
               end
            end

            xigma3_1=-xigma3;
            xigma3=xigma3_1;

            a7_3=[];
            a7_4=[];
            d7_3=[];
            d7_4=[];
            n_tr7_3=[];
            n_tr7_4=[];

            n_tr7_3=zeros(8,3);
            n_tr7_4=zeros(8,3);
            if n_tr5(9,3)==0&n_tr5(10,3)==0
               n_tr7_3(1:2,1:3)=zeros(2,3);
               n_tr7_4(1:2,1:3)=zeros(2,3);
            end
            if n_tr5(9,3)==0&n_tr5(10,3)>0
               n_tr7_3(1,1:3)=zeros(1,3);
               n_tr7_3(2,1)=1;
               n_tr7_3(2,2)=floor(n_tr5(10,3)/2);
               n_tr7_3(2,3)=n_tr7_3(2,2);
               a7_3=a5(n_tr5(10,1):n_tr5(10,1)+floor(n_tr5(10,3)/2)-1,:);
               d7_3(:,1)=-ones(floor(n_tr5(10,3)/2),1);

               n_tr7_4(1,1:3)=zeros(1,3);
               n_tr7_4(2,1)=1;
               n_tr7_4(2,2)=n_tr5(10,3)-floor(n_tr5(10,3)/2);
               n_tr7_4(2,3)=n_tr7_4(2,2);
               a7_4=a5(n_tr5(10,1)+floor(n_tr5(10,3)/2):n_tr5(10,2),:);
               d7_4(:,1)=-ones(n_tr5(10,3)-floor(n_tr5(10,3)/2),1);
            end
            if n_tr5(9,3)>0&n_tr5(10,3)==0
               if n_tr5(9,3)>1
                  n_tr7_3(1,1)=1;
                  n_tr7_3(1,2)=floor(n_tr5(9,3)/2);
                  n_tr7_3(1,3)=n_tr7_3(1,2);
                  a7_3=a5(n_tr5(9,1):n_tr5(9,1)+floor(n_tr5(9,3)/2)-1,:);
                  d7_3(:,1)=ones(floor(n_tr5(9,3)/2),1);
               else
                  n_tr7_3(1,1:3)=zeros(1,3);
               end
               n_tr7_3(2,1:3)=zeros(1,3);

               n_tr7_4(1,1)=1;
               n_tr7_4(1,2)=n_tr5(9,3)-floor(n_tr5(9,3)/2);
               n_tr7_4(1,3)=n_tr7_4(1,2);
               n_tr7_4(2,1:3)=zeros(1,3);
               a7_4=a5(n_tr5(9,1)+floor(n_tr5(9,3)/2):n_tr5(9,2),:);
               d7_4(:,1)=ones(n_tr5(9,3)-floor(n_tr5(9,3)/2),1);
            end
            if (n_tr5(9,3)>0&n_tr5(10,3)>0)&minnumerror6(5)==0
               if n_tr5(9,3)>1
                  n_tr7_3(1,1)=1;
                  n_tr7_3(1,2)=floor(n_tr5(9,3)/2);
                  n_tr7_3(1,3)=n_tr7_3(1,2);
                  a7_3=a5(n_tr5(9,1):n_tr5(9,1)+floor(n_tr5(9,3)/2)-1,:);
                  d7_3(:,1)=ones(floor(n_tr5(9,3)/2),1);
               else
                  n_tr7_3(1,1:3)=zeros(1,3);
               end
               n_tr7_3(2,1)=n_tr7_3(1,2)+1;
               n_tr7_3(2,2)=n_tr7_3(1,2)+floor(n_tr5(10,3)/2);
               n_tr7_3(2,3)=floor(n_tr5(10,3)/2);

               a7_3(n_tr7_3(2,1):n_tr7_3(2,2),:)=a5(n_tr5(10,1):n_tr5(10,1)+floor(n_tr5(10,3)/2)-1,:);
               d7_3(n_tr7_3(2,1):n_tr7_3(2,2),1)=-ones(floor(n_tr5(10,3)/2),1);

               n_tr7_4(1,1)=1;
               n_tr7_4(1,2)=n_tr5(9,3)-floor(n_tr5(9,3)/2);
               n_tr7_4(1,3)=n_tr7_4(1,2);
               n_tr7_4(2,1)=n_tr7_4(1,2)+1;
               n_tr7_4(2,2)=n_tr7_4(1,2)+n_tr5(10,3)-floor(n_tr5(10,3)/2);
               n_tr7_4(2,3)=n_tr5(10,3)-floor(n_tr5(10,3)/2);

               a7_4=a5(n_tr5(9,1)+floor(n_tr5(9,3)/2):n_tr5(9,2),:);
               a7_4(n_tr7_4(2,1):n_tr7_4(2,2),:)=a5(n_tr5(9,2)+floor(n_tr5(10,3)/2)+1:n_tr5(10,2),:);
               d7_4(:,1)=ones(n_tr5(9,3)-floor(n_tr5(9,3)/2),1);
               d7_4(n_tr7_4(2,1):n_tr7_4(2,2),1)=-ones(n_tr5(10,3)-floor(n_tr5(10,3)/2),1);
            end

            if (n_tr5(9,3)>0&n_tr5(10,3)>0)&minnumerror6(5)>0
               n_tr6=[];
               n_tr6(1,1)=1;
               n_tr6(1,2)=n_tr5(9,3);
               n_tr6(1,3)=n_tr5(9,3);
               n_tr6(1,4)=1;
               n_tr6(2,1)=n_tr5(9,3)+1;
               n_tr6(2,2)=n_tr5(9,3)+n_tr5(10,3);
               n_tr6(2,3)=n_tr5(10,3);
               n_tr6(2,4)=2;
               a00=[];
               a00=a5(n_tr5(9,1):n_tr5(10,2),:);

              [a71,d71,n_tr71]=sample_decomposition_8(a00,n_tr6,w5,theta5,xigma3,dim);
               n_tr7_3(1:2,1:3)=n_tr71(1:2,1:3);

               if n_tr71(1,3)==0&n_tr71(2,3)>0
                  a7_3=a71(n_tr71(2,1):n_tr71(2,2),:);
                  d7_3=d71(n_tr71(2,1):n_tr71(2,2),:);
               end
               if n_tr71(1,3)>0&n_tr71(2,3)==0
                  a7_3=a71(n_tr71(1,1):n_tr71(1,2),:);
                  d7_3=d71(n_tr71(1,1):n_tr71(1,2),:);
               end
               if n_tr71(1,3)>0&n_tr71(2,3)>0
                  a7_3=a71(n_tr71(1,1):n_tr71(2,2),:);
                  d7_3=d71(n_tr71(1,1):n_tr71(2,2),:);
               end

               if n_tr71(3,3)==0&n_tr71(4,3)==0
                  n_tr7_4(1:2,1:3)=zeros(2,3);
               end
               if n_tr71(3,3)==0&n_tr71(4,3)>0
                  n_tr7_4(1,1:3)=zeros(1,3);
                  n_tr7_4(2,1)=1;
                  n_tr7_4(2,3)=n_tr71(4,3);
                  n_tr7_4(2,2)=n_tr71(4,3);
                  a7_4=a71(n_tr71(4,1):n_tr71(4,2),:);
                  d7_4=d71(n_tr71(4,1):n_tr71(4,2),:);
               end
               if n_tr71(3,3)>0&n_tr71(4,3)==0
                  n_tr7_4(1,1)=1;
                  n_tr7_4(1,2)=n_tr71(3,3);
                  n_tr7_4(1,3)=n_tr71(3,3);
                  n_tr7_4(2,1:3)=zeros(1,3);
                  a7_4=a71(n_tr71(3,1):n_tr71(3,2),:);
                  d7_4=d71(n_tr71(3,1):n_tr71(3,2),:);
               end
               if n_tr71(3,3)>0&n_tr71(4,3)>0
                  n_tr7_4(1,1)=1;
                  n_tr7_4(1,2)=n_tr71(3,3);
                  n_tr7_4(2,1)=n_tr71(3,3)+1;
                  n_tr7_4(2,2)=n_tr71(3,3)+n_tr71(4,3);
                  n_tr7_4(1:2,3)=n_tr71(3:4,3);
                  a7_4=a71(n_tr71(3,1):n_tr71(4,2),:);
                  d7_4=d71(n_tr71(3,1):n_tr71(4,2),:);
               end
            end

            for ii=6:8,
                n_tr7_3(2*ii-9:2*ii-8,1:3)=zeros(2,3);
                n_tr7_4(2*ii-9:2*ii-8,1:3)=zeros(2,3);
                if n_tr5(2*ii-1,3)==0&n_tr5(2*ii,3)==0
                   n_tr7_3(2*ii-9:2*ii-8,1:3)=zeros(2,3);
                   n_tr7_4(2*ii-9:2*ii-8,1:3)=zeros(2,3);
                end
                if n_tr5(2*ii-1,3)==0&n_tr5(2*ii,3)>0
                   n_tr7_3(2*ii-9,1:3)=zeros(1,3);
                   n_tr7_3(2*ii-8,1)=sum(n_tr7_3(1:2*(ii-5),3))+1;
                   n_tr7_3(2*ii-8,2)=sum(n_tr7_3(1:2*(ii-5),3))+floor(n_tr5(2*ii,3)/2);
                   n_tr7_3(2*ii-8,3)=floor(n_tr5(2*ii,3)/2);
                   a7_3(n_tr7_3(2*ii-8,1):n_tr7_3(2*ii-8,2),:)=a5(sum(n_tr5(1:2*(ii-1),3))+1:sum(n_tr5(1:2*(ii-1),3))+floor(n_tr5(2*ii,3)/2),:);
                   d7_3(n_tr7_3(2*ii-8,1):n_tr7_3(2*ii-8,2),1)=-ones(floor(n_tr5(2*ii,3)/2),1);

                   n_tr7_4(2*ii-9,1:3)=zeros(1,3);
                   n_tr7_4(2*ii-8,1)=sum(n_tr7_4(1:2*(ii-5),3))+1;
                   n_tr7_4(2*ii-8,2)=sum(n_tr7_4(1:2*(ii-5),3))+n_tr5(2*ii,3)-floor(n_tr5(2*ii,3)/2);
                   n_tr7_4(2*ii-8,3)=n_tr5(2*ii,3)-floor(n_tr5(2*ii,3)/2);
                   a7_4(n_tr7_4(2*ii-8,1):n_tr7_4(2*ii-8,2),:)=a5(sum(n_tr5(1:2*(ii-1),3))+floor(n_tr5(2*ii,3)/2)+1:sum(n_tr5(1:2*(ii-1),3))+n_tr5(2*ii,3),:);
                   d7_4(n_tr7_4(2*ii-8,1):n_tr7_4(2*ii-8,2),1)=-ones(n_tr5(2*ii,3)-floor(n_tr5(2*ii,3)/2),1);
                end
                if n_tr5(2*ii-1,3)>0&n_tr5(2*ii,3)==0
                   n_tr7_3(2*ii-9,1)=sum(n_tr7_3(1:2*(ii-5),3))+1;
                   n_tr7_3(2*ii-9,2)=sum(n_tr7_3(1:2*(ii-5),3))+floor(n_tr5(2*ii-1,3)/2);
                   n_tr7_3(2*ii-9,3)=floor(n_tr5(2*ii-1,3)/2);
                   a7_3(n_tr7_3(2*ii-9,1):n_tr7_3(2*ii-9,2),:)=a5(sum(n_tr5(1:2*(ii-1),3))+1:sum(n_tr5(1:2*(ii-1),3))+floor(n_tr5(2*ii-1,3)/2),:);
                   d7_3(n_tr7_3(2*ii-9,1):n_tr7_3(2*ii-9,2),1)=ones(floor(n_tr5(2*ii-1,3)/2),1);
                   n_tr7_3(2*ii-8,1:3)=zeros(1,3);

                   n_tr7_4(2*ii-9,1)=sum(n_tr7_4(1:2*(ii-5),3))+1;
                   n_tr7_4(2*ii-9,2)=sum(n_tr7_4(1:2*(ii-5),3))+n_tr5(2*ii-1,3)-floor(n_tr5(2*ii-1,3)/2);
                   n_tr7_4(2*ii-9,3)=n_tr5(2*ii-1,3)-floor(n_tr5(2*ii-1,3)/2);
                   a7_4(n_tr7_4(2*ii-9,1):n_tr7_4(2*ii-9,2),:)=a5(sum(n_tr5(1:2*(ii-1),3))+floor(n_tr5(2*ii-1,3)/2)+1:sum(n_tr5(1:2*(ii-1),3))+n_tr5(2*ii-1,3),:);
                   d7_4(n_tr7_4(2*ii-9,1):n_tr7_4(2*ii-9,2),1)=ones(n_tr5(2*ii-1,3)-floor(n_tr5(2*ii-1,3)/2),1);
                   n_tr7_4(2*ii-8,1:3)=zeros(1,3);
                end
                if n_tr5(2*ii-1,3)>0&n_tr5(2*ii,3)>0&minnumerror6(ii)==0
                   n_tr7_3(2*ii-9,1)=sum(n_tr7_3(1:2*(ii-5),3))+1;
                   n_tr7_3(2*ii-9,2)=sum(n_tr7_3(1:2*(ii-5),3))+floor(n_tr5(2*ii-1,3)/2);
                   n_tr7_3(2*ii-9,3)=floor(n_tr5(2*ii-1,3)/2);
                   n_tr7_3(2*ii-8,1)=sum(n_tr7_3(1:2*ii-9,3))+1;
                   n_tr7_3(2*ii-8,2)=sum(n_tr7_3(1:2*ii-9,3))+floor(n_tr5(2*ii,3)/2);
                   n_tr7_3(2*ii-8,3)=floor(n_tr5(2*ii,3)/2);
                   a7_3(n_tr7_3(2*ii-9,1):n_tr7_3(2*ii-9,2),:)=a5(sum(n_tr5(1:2*(ii-1),3))+1:sum(n_tr5(1:2*(ii-1),3))+floor(n_tr5(2*ii-1,3)/2),:);
                   d7_3(n_tr7_3(2*ii-9,1):n_tr7_3(2*ii-9,2),1)=ones(floor(n_tr5(2*ii-1,3)/2),1);
                   a7_3(n_tr7_3(2*ii-8,1):n_tr7_3(2*ii-8,2),:)=a5(sum(n_tr5(1:2*ii-1,3))+1:sum(n_tr5(1:2*ii-1,3))+floor(n_tr5(2*ii,3)/2),:);
                   d7_3(n_tr7_3(2*ii-8,1):n_tr7_3(2*ii-8,2),1)=-ones(floor(n_tr5(2*ii,3)/2),1);

                   n_tr7_4(2*ii-9,1)=sum(n_tr7_4(1:2*(ii-5),3))+1;
                   n_tr7_4(2*ii-9,2)=sum(n_tr7_4(1:2*(ii-5),3))+n_tr5(2*ii-1,3)-floor(n_tr5(2*ii-1,3)/2);
                   n_tr7_4(2*ii-9,3)=n_tr5(2*ii-1,3)-floor(n_tr5(2*ii-1,3)/2);
                   n_tr7_4(2*ii-8,1)=sum(n_tr7_4(1:2*ii-9,3))+1;
                   n_tr7_4(2*ii-8,2)=sum(n_tr7_4(1:2*ii-9,3))+n_tr5(2*ii,3)-floor(n_tr5(2*ii,3)/2);
                   n_tr7_4(2*ii-8,3)=n_tr5(2*ii,3)-floor(n_tr5(2*ii,3)/2);
                   a7_4(n_tr7_4(2*ii-9,1):n_tr7_4(2*ii-9,2),:)=a5(sum(n_tr5(1:2*(ii-1),3))+floor(n_tr5(2*ii-1,3)/2)+1:n_tr5(2*ii-1,2),:);
                   d7_4(n_tr7_4(2*ii-9,1):n_tr7_4(2*ii-9,2),1)=ones(n_tr5(2*ii-1,3)-floor(n_tr5(2*ii-1,3)/2),1);
                   a7_4(n_tr7_4(2*ii-8,1):n_tr7_4(2*ii-8,2),:)=a5(sum(n_tr5(1:2*ii-1,3))+floor(n_tr5(2*ii,3)/2)+1:n_tr5(2*ii,2),:);
                   d7_4(n_tr7_4(2*ii-8,1):n_tr7_4(2*ii-8,2),1)=-ones(n_tr5(2*ii,3)-floor(n_tr5(2*ii,3)/2),1);
                end

                if n_tr5(2*ii-1,3)>0&n_tr5(2*ii,3)>0&minnumerror6(ii)>0
                   n_tr6=[];
                   n_tr6(1,1)=1;
                   n_tr6(1,2)=n_tr5(2*ii-1,3);
                   n_tr6(1,3)=n_tr5(2*ii-1,3);
                   n_tr6(1,4)=1;
                   n_tr6(2,1)=n_tr5(2*ii-1,3)+1;
                   n_tr6(2,2)=n_tr5(2*ii-1,3)+n_tr5(2*ii,3);
                   n_tr6(2,3)=n_tr5(2*ii,3);
                   n_tr6(2,4)=2;
                   a00=[];
                   a00=a5(n_tr5(2*ii-1,1):n_tr5(2*ii,2),:);
                  [a71,d71,n_tr71]=sample_decomposition_8(a00,n_tr6,w5,theta5,xigma3,dim);

                   n_tr7_3(2*ii-9:2*ii-8,1:2)=sum(n_tr7_3(1:2*(ii-5),3))+n_tr71(1:2,1:2);
                   n_tr7_3(2*ii-9:2*ii-8,3)=n_tr71(1:2,3);
                   if n_tr71(1,3)==0&n_tr71(2,3)>0
                      a7_3(n_tr7_3(2*ii-8,1):n_tr7_3(2*ii-8,2),:)=a71(n_tr71(2,1):n_tr71(2,2),:);
                      d7_3(n_tr7_3(2*ii-8,1):n_tr7_3(2*ii-8,2),1)=d71(n_tr71(2,1):n_tr71(2,2),:);
                   end
                   if n_tr71(1,3)>0&n_tr71(2,3)==0
                      a7_3(n_tr7_3(2*ii-9,1):n_tr7_3(2*ii-9,2),:)=a71(n_tr71(1,1):n_tr71(1,2),:);
                      d7_3(n_tr7_3(2*ii-9,1):n_tr7_3(2*ii-9,2),1)=d71(n_tr71(1,1):n_tr71(1,2),:);
                   end
                   if n_tr71(1,3)>0&n_tr71(2,3)>0
                      a7_3(n_tr7_3(2*ii-9,1):n_tr7_3(2*ii-8,2),:)=a71(n_tr71(1,1):n_tr71(2,2),:);
                      d7_3(n_tr7_3(2*ii-9,1):n_tr7_3(2*ii-8,2),1)=d71(n_tr71(1,1):n_tr71(2,2),:);
                   end

                   if n_tr71(3,3)==0&n_tr71(4,3)==0
                      n_tr7_4(2*ii-9:2*ii-8,1:3)=zeros(2,3);
                   end
                   if n_tr71(3,3)==0&n_tr71(4,3)>0
                      n_tr7_4(2*ii-9,1:3)=zeros(1,3);
                      n_tr7_4(2*ii-8,3)=n_tr71(4,3);
                      n_tr7_4(2*ii-8,1)=sum(n_tr7_4(1:2*ii-9,3))+1;
                      n_tr7_4(2*ii-8,2)=sum(n_tr7_4(1:2*ii-9,3))+n_tr71(4,3);
                      a7_4(n_tr7_4(2*ii-8,1):n_tr7_4(2*ii-8,2),:)=a71(n_tr71(4,1):n_tr71(4,2),:);
                      d7_4(n_tr7_4(2*ii-8,1):n_tr7_4(2*ii-8,2),1)=d71(n_tr71(4,1):n_tr71(4,2),:);
                   end
                   if n_tr71(3,3)>0&n_tr71(4,3)==0
                      n_tr7_4(2*ii-9,1)=sum(n_tr7_4(1:2*(ii-5),3))+1;
                      n_tr7_4(2*ii-9,2)=sum(n_tr7_4(1:2*(ii-5),3))+n_tr71(3,3);
                      n_tr7_4(2*ii-9,3)=n_tr71(3,3);
                      n_tr7_4(2*ii-8,1:3)=zeros(1,3);
                      a7_4(n_tr7_4(2*ii-9,1):n_tr7_4(2*ii-9,2),:)=a71(n_tr71(3,1):n_tr71(3,2),:);
                      d7_4(n_tr7_4(2*ii-9,1):n_tr7_4(2*ii-9,2),1)=d71(n_tr71(3,1):n_tr71(3,2),:);
                   end
                   if n_tr71(3,3)>0&n_tr71(4,3)>0
                      n_tr7_4(2*ii-9:2*ii-8,3)=n_tr71(3:4,3);
                      n_tr7_4(2*ii-9,1)=sum(n_tr7_4(1:2*ii-10,3))+1;
                      n_tr7_4(2*ii-9,2)=sum(n_tr7_4(1:2*ii-10,3))+n_tr71(3,3);
                      n_tr7_4(2*ii-8,1)=sum(n_tr7_4(1:2*ii-9,3))+1;
                      n_tr7_4(2*ii-8,2)=sum(n_tr7_4(1:2*ii-9,3))+n_tr71(4,3);
                      a7_4(n_tr7_4(2*ii-9,1):n_tr7_4(2*ii-8,2),:)=a71(n_tr71(3,1):n_tr71(4,2),:);
                      d7_4(n_tr7_4(2*ii-9,1):n_tr7_4(2*ii-8,2),1)=d71(n_tr71(3,1):n_tr71(4,2),:);
                   end
               end
            end

            n_tr7=n_tr7_1;
            n_tr7(9:16,:)=n_tr7_2;
            n_tr7(17:24,:)=n_tr7_3;
            n_tr7(25:32,:)=n_tr7_4;
            n_tr7(9:16,1:2)=n_tr7_2(:,1:2)+sum(n_tr7(1:8,3));
            n_tr7(17:24,1:2)=n_tr7_3(:,1:2)+sum(n_tr7(1:16,3));
            n_tr7(25:32,1:2)=n_tr7_4(:,1:2)+sum(n_tr7(1:24,3));

            for i=1:32,
                n_tr7(i,4)=i;
            end
            for i=1:32,
                if n_tr7(i,3)==0
                   n_tr7(i,1:2)=[0 0];
                end
            end

            a7=a7_1;
            d7=d7_1;
            for i=9:16,
                if n_tr7(i,1)>0
                   a7(n_tr7(i,1):sum(n_tr7(1:16,3)),:)=a7_2;
                   d7(n_tr7(i,1):sum(n_tr7(1:16,3)),:)=d7_2;
                   break;
                end
            end
            for i=17:24,
                if n_tr7(i,1)>0
                   a7(n_tr7(i,1):sum(n_tr7(1:24,3)),:)=a7_3;
                   d7(n_tr7(i,1):sum(n_tr7(1:24,3)),:)=d7_3;
                   break;
                end
            end
            for i=25:32,
                if n_tr7(i,1)>0
                   a7(n_tr7(i,1):sum(n_tr7(:,3)),:)=a7_4;
                   d7(n_tr7(i,1):sum(n_tr7(:,3)),:)=d7_4;
                   break;
                end
            end

            n_tr7
            fprintf(fid,'\n n_tr7=\n',n_tr7');
            fprintf(fid,'%g  %g  %g  %g\n',n_tr7');
            fprintf(fid,'\n');

%/////////////////////////////////////////////////////////////////////////////////////////////////
            minnumerror7_1=zeros(1,16);
            minnumerror7_2=zeros(16,2);
%/////////////////////////////////////////////////////////////////////////////////////////////////
            for i=1:4,
                if n_tr7(8*i-7,3)==0&n_tr7(8*i-6,3)>0
                   w7_1(4*i-3,:)=-w2(k0,:);
                   theta7_1(4*i-3)=-theta2(k0);
                   n0000(4*i-3)=0;
                end
                if n_tr7(8*i-7,3)>=0&n_tr7(8*i-6,3)==0
                   w7_1(4*i-3,:)=w2(k0,:);
                   theta7_1(4*i-3)=theta2(k0);
                   n0000(4*i-3)=0;
                end

                if n_tr7(8*i-5,3)==0&n_tr7(8*i-4,3)>=0
                   w7_1(4*i-2,:)=w2(k0,:);
                   theta7_1(4*i-2)=theta2(k0);
                   n0000(4*i-2)=0;
                end
                if n_tr7(8*i-5,3)>0&n_tr7(8*i-4,3)==0
                   w7_1(4*i-2,:)=-w2(k0,:);
                   theta7_1(4*i-2)=-theta2(k0);
                   n0000(4*i-2)=0;
                end

                if n_tr7(8*i-3,3)==0&n_tr7(8*i-2,3)>0
                   w7_1(4*i-1,:)=-w2(k0,:);
                   theta7_1(4*i-1)=mu01*w2(k0,:)';
                   n0000(4*i-1)=0;
                end
                if n_tr7(8*i-3,3)>=0&n_tr7(8*i-2,3)==0
                   w7_1(4*i-1,:)=w2(k0,:);
                   theta7_1(4*i-1)=-mu01*w2(k0,:)';
                   n0000(4*i-1)=0;
                end

                if n_tr7(8*i-1,3)==0&n_tr7(8*i,3)>=0
                   w7_1(4*i,:)=w2(k0,:);
                   theta7_1(4*i)=-mu02*w2(k0,:)';
                   n0000(4*i)=0;
                end
                if n_tr7(8*i-1,3)>0&n_tr7(8*i,3)==0
                   w7_1(4*i,:)=-w2(k0,:);
                   theta7_1(4*i)=mu02*w2(k0,:)';
                   n0000(4*i)=0;
                end
            end

            minnumerror66=minnumerror6(1:4);
            minnumerror66(5:8)=minnumerror6(1:4);
            minnumerror66(9:12)=minnumerror6(5:8);
            minnumerror66(13:16)=minnumerror6(5:8);
            ww6=w66(1:4,:);
            ww6(5:8,:)=w66(1:4,:);
            ww6(9:12,:)=w66(5:8,:);
            ww6(13:16,:)=w66(5:8,:);
            thethe6=theta66(1:4);
            thethe6(5:8)=theta66(1:4);
            thethe6(9:12)=theta66(5:8);
            thethe6(13:16)=theta66(5:8);

            for i=1:16,
                p2=n_tr7(2*i-1,3);
                q2=n_tr7(2*i,3);
                a22=[];
                d22=[];
                if p2>0&q2>0&minnumerror66(i)>0
                   a22=a7(n_tr7(2*i-1,1):n_tr7(2*i,2),:);
                   d22=d7(n_tr7(2*i-1,1):n_tr7(2*i,2),:);
                  [w3,theta3,minnumerror3,n]=calculating_weigh_theta_6(a22,dim,d22,p0,q0,p2,q2,fid);
                   fprintf('num(%d,%d,%d)=%d\n\n',r,t,i,n);
                   fprintf(fid,'num(%d,%d,%d)=%d\n\n',r,t,i,n);
                   w7_1(i,:)=w3';
                   theta7_1(i)=theta3;
                   minnumerror7_1(i)=sum(minnumerror3(2:3));
                   minnumerror7_2(i,:)=minnumerror3(2:3);
                   n0000(i)=n;
                end
                if p2>0&q2>0&minnumerror66(i)==0
                   w7_1(i,:)=ww6(i,:);
                   theta7_1(i)=thethe6(i);
                   minnumerror7_1(i)=0;
                   minnumerror7_2(i,:)=[0,0];
                   n0000(i)=0;
                end
            end
            n000(k0)=n000(k0)+sum(n0000);

            numtr0=zeros(1,2);
            numtr1=zeros(16,2);
            for i=1:16,
                p2=n_tr7(2*i-1,3);
                q2=n_tr7(2*i,3);
                a22=[];
                d22=[];
                re2=[];
                if p2>0&q2>0
                   a22=a7(n_tr7(2*i-1,1):n_tr7(2*i,2),:);
                   d22=d7(n_tr7(2*i-1,1):n_tr7(2*i,2),:);
                end
                if p2==0&q2>0
                   a22=a7(n_tr7(2*i,1):n_tr7(2*i,2),:);
                   d22=d7(n_tr7(2*i,1):n_tr7(2*i,2),:);
                end
                if p2>0&q2==0
                   a22=a7(n_tr7(2*i-1,1):n_tr7(2*i-1,2),:);
                   d22=d7(n_tr7(2*i-1,1):n_tr7(2*i-1,2),:);
                end
                if p2>0|q2>0
                   re2=a22(:,2:dim+1)*(w7_1(i,:))'+ones(p2+q2,1)*theta7_1(i);
                end
                if p2>0
                   for ii=1:p2,
                       if re2(ii)<0
                          numtr0(1)=numtr0(1)+1;
                          numtr1(i,1)=numtr1(i,1)+1;
                       end
                   end
                end
                if q2>0
                   for ii=p2+1:p2+q2,
                       if re2(ii)>0
                          numtr0(2)=numtr0(2)+1;
                          numtr1(i,2)=numtr1(i,2)+1;
                       end
                   end
                end
            end
            numtr0
            numtr1'
            sum(sum(numtr1))
            fprintf(fid,'\n number of error samples is\n');
            fprintf(fid,'%g  %g\n',numtr0');
            fprintf(fid,'\n%g  %g\n',numtr1');
            fprintf(fid,'\n');

            minnumerror7=minnumerror7_1;
            fprintf('minnumerror7=%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, sum=%g\n',minnumerror7,sum(minnumerror7));
            fprintf(fid,'minnumerror7=%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, sum=%g\n',minnumerror7,sum(minnumerror7));

            s2=sum(sign0)-11;
            w8(s2,:)=w2(k0,:);
            theta8(s2)=theta2(k0);

            sign0(k0)=22;
            minnumerror8_1=[];

            w8_1=w7_1;
            theta8_1=theta7_1;
            minnumerror8_1=minnumerror7_1;

            if minnumerror7_2(1,1)+minnumerror7_2(1,2)>=n_tr7(2,3)
               w8_1(1,:)=w66(1,:);
               theta8_1(1)=theta66(1);
               minnumerror8_1(1)=n_tr7(2,3);
            end
            if minnumerror7_2(5,1)+minnumerror7_2(5,2)>=n_tr7(10,3)
               w8_1(5,:)=w66(1,:);
               theta8_1(5)=theta66(1);
               minnumerror8_1(5)=n_tr7(10,3);
            end
            if minnumerror6(1)==0
               w8_1(1,:)=w66(1,:);
               theta8_1(1)=theta66(1);
               w8_1(5,:)=w66(1,:);
               theta8_1(5)=theta66(1);
               minnumerror8_1(1)=0;
               minnumerror8_1(5)=0;
            end
            if minnumerror7_2(9,1)+minnumerror7_2(9,2)>=n_tr7(18,3)
               w8_1(9,:)=w66(1,:);
               theta8_1(9)=theta66(1);
               minnumerror8_1(9)=n_tr7(18,3);
            end
            if minnumerror7_2(13,1)+minnumerror7_2(13,2)>=n_tr7(26,3)
               w8_1(13,:)=w66(1,:);
               theta8_1(13)=theta66(1);
               minnumerror8_1(13)=n_tr7(26,3);
            end
            if minnumerror6(5)==0
               w8_1(9,:)=w66(5,:);
               theta8_1(9)=theta66(5);
               w8_1(13,:)=w66(5,:);
               theta8_1(13)=theta66(5);
               minnumerror8_1(9)=0;
               minnumerror8_1(13)=0;
            end

            if minnumerror7_2(2,1)+minnumerror7_2(2,2)>=n_tr7(3,3)
               w8_1(2,:)=w66(2,:);
               theta8_1(2)=theta66(2);
               minnumerror8_1(2)=n_tr7(3,3);
            end
            if minnumerror7_2(6,1)+minnumerror7_2(6,2)>=n_tr7(11,3)
               w8_1(6,:)=w66(2,:);
               theta8_1(6)=theta66(2);
               minnumerror8_1(6)=n_tr7(11,3);
            end
            if minnumerror6(2)==0
               w8_1(2,:)=w66(2,:);
               theta8_1(2)=theta66(2);
               w8_1(6,:)=w66(2,:);
               theta8_1(6)=theta66(2);
               minnumerror8_1(2)=0;
               minnumerror8_1(6)=0;
            end
            if minnumerror7_2(10,1)+minnumerror7_2(10,2)>=n_tr7(19,3)
               w8_1(10,:)=w66(6,:);
               theta8_1(10)=theta66(6);
               minnumerror8_1(10)=n_tr7(19,3);
            end
            if minnumerror7_2(14,1)+minnumerror7_2(14,2)>=n_tr7(27,3)
               w8_1(14,:)=w66(6,:);
               theta8_1(14)=theta66(6);
               minnumerror8_1(14)=n_tr7(27,3);
            end
            if minnumerror6(6)==0
               w8_1(10,:)=w66(6,:);
               theta8_1(10)=theta66(6);
               w8_1(14,:)=w66(6,:);
               theta8_1(14)=theta66(6);
               minnumerror8_1(10)=0;
               minnumerror8_1(14)=0;
            end

            if minnumerror7_2(3,1)+minnumerror7_2(3,2)>=n_tr7(6,3)
               w8_1(3,:)=w66(3,:);
               theta8_1(3)=theta66(3);
               minnumerror8_1(3)=n_tr7(6,3);
            end
            if minnumerror7_2(7,1)+minnumerror7_2(7,2)>=n_tr7(14,3)
               w8_1(7,:)=w66(3,:);
               theta8_1(7)=theta66(3);
               minnumerror8_1(7)=n_tr7(14,3);
            end
            if minnumerror6(3)==0
               w8_1(3,:)=w66(3,:);
               theta8_1(3)=theta66(3);
               w8_1(7,:)=w66(3,:);
               theta8_1(7)=theta66(3);
               minnumerror8_1(3)=0;
               minnumerror8_1(7)=0;
            end
            if minnumerror7_2(11,1)+minnumerror7_2(11,2)>=n_tr7(22,3)
               w8_1(11,:)=w66(7,:);
               theta8_1(11)=theta66(7);
               minnumerror8_1(11)=n_tr7(22,3);
            end
            if minnumerror7_2(15,1)+minnumerror7_2(15,2)>=n_tr7(30,3)
               w8_1(15,:)=w66(7,:);
               theta8_1(15)=theta66(7);
               minnumerror8_1(15)=n_tr7(30,3);
            end
            if minnumerror6(7)==0
               w8_1(11,:)=w66(7,:);
               theta8_1(11)=theta66(7);
               w8_1(15,:)=w66(7,:);
               theta8_1(15)=theta66(7);
               minnumerror8_1(11)=0;
               minnumerror8_1(15)=0;
            end

            if minnumerror7_2(4,1)+minnumerror7_2(4,2)>=n_tr7(7,3)
               w8_1(4,:)=w66(4,:);
               theta8_1(4)=theta66(4);
               minnumerror8_1(4)=n_tr7(7,3);
            end
            if minnumerror7_2(8,1)+minnumerror7_2(8,2)>=n_tr7(15,3)
               w8_1(8,:)=w66(4,:);
               theta8_1(8)=theta66(4);
               minnumerror8_1(8)=n_tr7(15,3);
            end
            if minnumerror6(4)==0
               w8_1(4,:)=w66(4,:);
               theta8_1(4)=theta66(4);
               w8_1(8,:)=w66(4,:);
               theta8_1(8)=theta66(4);
               minnumerror8_1(4)=0;
               minnumerror8_1(8)=0;
            end
            if minnumerror7_2(12,1)+minnumerror7_2(12,2)>=n_tr7(23,3)
               w8_1(12,:)=w66(8,:);
               theta8_1(12)=theta66(8);
               minnumerror8_1(12)=n_tr7(23,3);
            end
            if minnumerror7_2(16,1)+minnumerror7_2(16,2)>=n_tr7(31,3)
               w8_1(16,:)=w66(8,:);
               theta8_1(16)=theta66(8);
               minnumerror8_1(16)=n_tr7(31,3);
            end
            if minnumerror6(8)==0
               w8_1(12,:)=w66(8,:);
               theta8_1(12)=theta66(8);
               w8_1(16,:)=w66(8,:);
               theta8_1(16)=theta66(8);
               minnumerror8_1(12)=0;
               minnumerror8_1(16)=0;
            end

%/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

            w8(s2+1,:)=w2(k0,:);
            theta8(s2+1)=-mu01*w2(k0,:)';
            w8(s2+2,:)=w2(k0,:);
            theta8(s2+2)=-mu02*w2(k0,:)';
            w8(s2+3,:)=w5;
            theta8(s2+3)=theta5;
            w8(s2+4,:)=w5;
            theta8(s2+4)=theta5+xigma3;
            w8(s2+5,:)=w5;
            theta8(s2+5)=theta5-xigma3;
            w8(s2+6:s2+21,:)=w8_1;
            theta8(s2+6:s2+21)=theta8_1;

            minnumerror8=minnumerror8_1;
            fprintf('minnumerror8=%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n',minnumerror8');
            fprintf('sum(minnumerror8)=%g\n',sum(minnumerror8));
            fprintf(fid,'minnumerror8=%g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g, %g\n',minnumerror8');
            fprintf(fid,'sum(minnumerror8)=%g\n',sum(minnumerror8));

%//////////////////////////////////////////////////
            numtr0=zeros(1,2);
            numtr1=zeros(16,2);
            kk=0;
            for i=1:16,
                p2=n_tr7(2*i-1,3);
                q2=n_tr7(2*i,3);
                a22=[];
                d22=[];
                re2=[];
                if p2>0&q2>0
                   a22=a7(n_tr7(2*i-1,1):n_tr7(2*i,2),:);
                   d22=d7(n_tr7(2*i-1,1):n_tr7(2*i,2),:);
                end
                if p2==0&q2>0
                   a22=a7(n_tr7(2*i,1):n_tr7(2*i,2),:);
                   d22=d7(n_tr7(2*i,1):n_tr7(2*i,2),:);
                end
                if p2>0&q2==0
                   a22=a7(n_tr7(2*i-1,1):n_tr7(2*i-1,2),:);
                   d22=d7(n_tr7(2*i-1,1):n_tr7(2*i-1,2),:);
                end
                if p2>0|q2>0
                   re2=a22(:,2:dim+1)*(w8_1(i,:))'+ones(p2+q2,1)*theta8_1(i);
                end
                if p2>0
                   for ii=1:p2,
                       if re2(ii)<0
                          numtr0(1)=numtr0(1)+1;
                          numtr1(i,1)=numtr1(i,1)+1;
                       end
                   end
                end
                if q2>0
                   for ii=p2+1:p2+q2,
                       if re2(ii)>0
                          numtr0(2)=numtr0(2)+1;
                          numtr1(i,2)=numtr1(i,2)+1;
                       end
                   end
                end
            end
            numtr0
            numtr1'
            sum(sum(numtr1))
            minnumerror2(k0)=sum(numtr0);
            fprintf(fid,'\nNumber of error samples is\n');
            fprintf(fid,'%g  %g\n\n',numtr0');
            fprintf(fid,'%g  %g\n',numtr1');
            fprintf(fid,'\n');
        end
        numtr_0(4)=numtr0(1);
        tacc_index(4)=sum(numtr0);

%%////////////////////////////////////////////////////////////////////////////
        pe=numtr0(1)/p0;
        qe=numtr0(2)/q0;
        if ir<=19
           aver_index(4)=(numtr0(1)+numtr0(2))/(p0+q0);
        elseif ir<=49
           aver_index(4)=(pe+qe)/2;
        elseif ir<=99.00
           if numtr0(1)>0&numtr0(2)>0
              aver_index(4)=(1-sqrt((1-pe)*(1-qe)));
           elseif numtr0(1)>0&numtr0(2)==0
              aver_index(4)=(1-sqrt(0.0001*(1-pe)));
           elseif numtr0(1)==0&numtr0(2)>0
              aver_index(4)=(1-sqrt(0.0001*(1-qe)));
           else
              aver_index(4)=0;
           end
        else
           if numtr0(1)>0&numtr0(2)>0
              aver_index(4)=(1/(1+alpha))*((pe+beta*qe-2*beta*pe*qe)/(2-pe-beta*qe)+alpha*pe);
           elseif numtr0(1)>0&numtr0(2)==0
              aver_index(4)=(1/(1+alpha))*((pe+beta*0.0001-2*beta*pe*0.0001)/(2-pe-beta*0.0001)+alpha*pe);
           elseif numtr0(1)==0&numtr0(2)>0
              aver_index(4)=(1/(1+alpha))*((0.0001+beta*qe-2*beta*0.0001*qe)/(2-0.0001-beta*qe)+alpha*pe);
           else
              aver_index(4)=0;
           end
        end
%%/////////////////////////////////////////////////////////////////////////////////////////////
aver_index
       [aver_min, imin]=min(aver_index);
       [tacc_min, imin1]=min(tacc_index);
        numtr_0
       [num_min,imin2]=min(numtr_0);
        imin
        imin1
        if imin<imin1
%           imin=imin1;
        end
        if imin>imin2&numtr_0(imin)>numtr_0(imin2)
%           imin=imin2;
        end
        aver_min
        imin

        if imin==1&aver_index(1)>0
           sign0(k0)=1;
           s3=sum(sign0);
           w4(s3,:)=zeros(1,dim);
           w4=w4(1:s3,:);
           theta4(s3)=0;
           theta4=theta4(1:s3);
           w6(s3,:)=zeros(1,dim);
           w6=w6(1:s3,:);
           theta6(s3)=0;
           theta6=theta6(1:s3);
           w8(s3,:)=zeros(1,dim);
           w8=w8(1:s3,:);
           theta8(s3)=0;
           theta8=theta8(1:s3);

           numtr0=[0 0];
           result02=a0(:,2:dim+1)*w2(k0,:)'+ones(p0+q0,1)*theta2(k0);
           re_mu(1)=mu01*w2(k0,:)'+theta2(k0);
           re_mu(2)=mu02*w2(k0,:)'+theta2(k0);
           numtr1=zeros(4,2);
           for i=1:p0,
               if result02(i)<0&result02(i)>=re_mu(2)
                  numtr0(1)=numtr0(1)+1;
                  numtr1(2,1)=numtr1(2,1)+1;
               end
               if result02(i)<0&result02(i)<re_mu(2)
                  numtr0(1)=numtr0(1)+1;
                  numtr1(4,1)=numtr1(4,1)+1;
               end
           end
           for i=p0+1:p0+q0,
               if result02(i)>=0&result02(i)<=re_mu(1)
                  numtr0(2)=numtr0(2)+1;
                  numtr1(1,2)=numtr1(1,2)+1;
               end
               if result02(i)>=0&result02(i)>re_mu(1)
                  numtr0(2)=numtr0(2)+1;
                  numtr1(3,2)=numtr1(3,2)+1;
               end
           end
           numtr0
           numtr1'
           sum(sum(numtr1))
        end

        if imin==2&aver_index(2)>0
           sign0(k0)=7;
           s3=sum(sign0)-6;
           w6(s3,:)=zeros(1,dim);
           w6=w6(1:s3,:);
           theta6(s3)=0;
           theta6=theta6(1:s3);
           w8(s3,:)=zeros(1,dim);
           w8=w8(1:s3,:);
           theta8(s3)=0;
           theta8=theta8(1:s3);

           numtr0=[0 0];
           re=a0(:,2:dim+1)*(w4(s3:s3+6,:))'+ones(p0+q0,1)*theta4(s3:s3+6);
           numtr1=zeros(4,2);
           for i=1:p0,
               if (re(i,1)>0)&(re(i,2)<0)&(re(i,4)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(1,1)=numtr1(1,1)+1;
               end
               if (re(i,1)<0)&(re(i,3)>0)&(re(i,5)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(2,1)=numtr1(2,1)+1;
               end
               if (re(i,1)>0)&(re(i,2)>0)&(re(i,6)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(3,1)=numtr1(3,1)+1;
               end
               if (re(i,1)<0)&(re(i,3)<0)&(re(i,7)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(4,1)=numtr1(4,1)+1;
               end
           end
           for i=p0+1:p0+q0,
               if (re(i,1)>0)&(re(i,2)<0)&(re(i,4)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(1,2)=numtr1(1,2)+1;
               end
               if (re(i,1)<0)&(re(i,3)>0)&(re(i,5)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(2,2)=numtr1(2,2)+1;
               end
               if (re(i,1)>0)&(re(i,2)>0)&(re(i,6)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(3,2)=numtr1(3,2)+1;
               end
               if (re(i,1)<0)&(re(i,3)<0)&(re(i,7)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(4,2)=numtr1(4,2)+1;
               end
           end
           numtr0
           numtr1'
           sum(sum(numtr1))
        end
        
        if imin==3&aver_index(3)>0
           sign0(k0)=12;
           s3=sum(sign0)-11;
           w8(s3,:)=zeros(1,dim);
           w8=w8(1:s3,:);
           theta8(s3)=0;
           theta8=theta8(1:s3);

           numtr0=[0 0];
           re1=a0(:,2:dim+1)*(w6(s3:s3+11,:))'+ones(p0+q0,1)*theta6(s3:s3+11);
           numtr1=zeros(8,2);
           for i=1:p0,
               if (re1(i,4)>0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,5)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(1,1)=numtr1(1,1)+1;
               end
               if (re1(i,4)>0)&(re1(i,1)<0)&(re1(i,3)>0)&(re1(i,6)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(2,1)=numtr1(2,1)+1;
               end
               if (re1(i,4)>0)&(re1(i,1)>0)&(re1(i,2)>0)&(re1(i,7)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(3,1)=numtr1(3,1)+1;
               end
               if (re1(i,4)>0)&(re1(i,1)<0)&(re1(i,3)<0)&(re1(i,8)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(4,1)=numtr1(4,1)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,9)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(5,1)=numtr1(5,1)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)<0)&(re1(i,3)>0)&(re1(i,10)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(6,1)=numtr1(6,1)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)>0)&(re1(i,2)>0)&(re1(i,11)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(7,1)=numtr1(7,1)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)<0)&(re1(i,3)<0)&(re1(i,12)<0)
                    numtr0(1)=numtr0(1)+1;
                   numtr1(8,1)=numtr1(8,1)+1;
               end
           end

           for i=p0+1:p0+q0,
               if (re1(i,4)>0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,5)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(1,2)=numtr1(1,2)+1;
               end
               if (re1(i,4)>0)&(re1(i,1)<0)&(re1(i,3)>0)&(re1(i,6)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(2,2)=numtr1(2,2)+1;
               end
               if (re1(i,4)>0)&(re1(i,1)>0)&(re1(i,2)>0)&(re1(i,7)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(3,2)=numtr1(3,2)+1;
               end
               if (re1(i,4)>0)&(re1(i,1)<0)&(re1(i,3)<0)&(re1(i,8)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(4,2)=numtr1(4,2)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,9)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(5,2)=numtr1(5,2)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)<0)&(re1(i,3)>0)&(re1(i,10)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(6,2)=numtr1(6,2)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)>0)&(re1(i,2)>0)&(re1(i,11)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(7,2)=numtr1(7,2)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)<0)&(re1(i,3)<0)&(re1(i,12)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(8,2)=numtr1(8,2)+1;
               end
           end
           numtr0
           numtr1'
           sum(sum(numtr1))
        end
        minnumerror2(k0)=sum(numtr0);

%//////////////////////////////////////////////////////////////
        fprintf('minnumerror2(%d,%d)=minnumerror2(%d)=%g\n',r,t,k0,minnumerror2(k0));
        fprintf(fid,'\nminnumerror2(%d,%d)=minnumerror2(%d)=%g\n',r,t,k0,minnumerror2(k0));
        fprintf(fid,'%g  ',minnumerror2);
        fprintf(fid,'\n');
        numtr00(k0)=sum(numtr0);
        floor(minnumerror2)
        t1=etime(clock,t0)
        fprintf(fid,'runtime(%d,%d)=%g\n\n',r,t,t1);
        total_time(k0)=t1;
        k0=k0+1;
    end
    if c==2
       break;
    end
end

t1=etime(clock,t0)
minnumerror2
n000
sum(n000)
total_time

fprintf(fid,'minnumerror2= ');
fprintf(fid,'%g  ',minnumerror2);
fprintf(fid,'\n');
fprintf(fid,'Total runtime=');
fprintf(fid,'%g  ',total_time);
fprintf(fid,'\n');
fprintf(fid,'Trial-and-error epochs=');
fprintf(fid,'%g  ',n000);
fprintf(fid,'\n');
fprintf(fid,'sum(n000)=%d\n',sum(n000));

fclose(fid);

if max(sign0)==1
   w4=zeros(1,dim);
   theta4=0;
end
if max(sign0)<=7
   w6=zeros(1,dim);
   theta6=0;
end

if max(sign0)<=12
   w8=zeros(1,dim);
   theta8=0;
end

fid=fopen('D:\highimbalance\letreca\inv_dc1\result\w0.txt','w');
    fprintf(fid,'%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n',w2');
fclose(fid);
fid=fopen('D:\highimbalance\letreca\inv_dc1\result\theta0.txt','w');
    fprintf(fid,'%.8f\n',theta2');
fclose(fid);
fid=fopen('D:\highimbalance\letreca\inv_dc1\result\w.txt','w');
    fprintf(fid,'%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n',w4');
fclose(fid);
fid=fopen('D:\highimbalance\letreca\inv_dc1\result\theta.txt','w');
    fprintf(fid,'%.8f\n',theta4');
fclose(fid);
fid=fopen('D:\highimbalance\letreca\inv_dc1\result\w1.txt','w');
    fprintf(fid,'%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n',w6');
fclose(fid);
fid=fopen('D:\highimbalance\letreca\inv_dc1\result\theta1.txt','w');
    fprintf(fid,'%.8f\n',theta6');
fclose(fid);
fid=fopen('D:\highimbalance\letreca\inv_dc1\result\w2.txt','w');
    fprintf(fid,'%.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f %.8f\n',w8');
fclose(fid);
fid=fopen('D:\highimbalance\letreca\inv_dc1\result\theta2.txt','w');
    fprintf(fid,'%.8f\n',theta8');
fclose(fid);
fid=fopen('D:\highimbalance\letreca\inv_dc1\result\sign0.txt','w');
    fprintf(fid,'%g\n',sign0');
fclose(fid);


%clear all;
%t0=clock;
% Total, average, geometric accuracy; AUC:
% For the rand_train1 set.
load D:\highimbalance\letreca\data_5rand\rand_train1.txt;
load D:\highimbalance\letreca\ntr.txt;
load D:\highimbalance\letreca\inv_dc1\result\w0.txt;
load D:\highimbalance\letreca\inv_dc1\result\w.txt;
load D:\highimbalance\letreca\inv_dc1\result\w1.txt;
load D:\highimbalance\letreca\inv_dc1\result\w2.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta0.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta1.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta2.txt;
load D:\highimbalance\letreca\inv_dc1\result\sign0.txt;

c=size(ntr,1);
dim=size(rand_train1,2)-2;
%dim=8;
k=0;

for r=1:c-1,
    for t=r+1:c,
        k=k+1;
        p0=ntr(r,3);
        q0=ntr(t,3);
        a=rand_train1(ntr(r,1):ntr(r,2),1:dim);
        a(p0+1:p0+q0,:)=rand_train1(ntr(t,1):ntr(t,2),1:dim);
        re0=a*(w0(k,:))'+ones(p0+q0,1)*theta0(k);
        kk=sum(sign0(1:k));
        if sign0(k)==7
           re=a*(w(kk-6:kk,:))'+ones(p0+q0,1)*theta(kk-6:kk)';
        end
        if sign0(k)==12
           re1=a*(w1(kk-11:kk,:))'+ones(p0+q0,1)*theta1(kk-11:kk)';
        end
        if sign0(k)==22
           re2=a*(w2(kk-21:kk,:))'+ones(p0+q0,1)*theta2(kk-21:kk)';
        end

% AUC
        y=[];
        [y(:,1),y(:,2)]=sort(re0,'descend');

        b=[0 0];
        j10=0;
        j00=0;
        j1=0;
        j0=0;

        for i=1:p0+q0,
            if y(i,2)<=p0
               j10=j10+1;
               j1=j10/p0;
            else
               j00=j00+1;
               j0=j00/q0;
            end
            b(i+1,1)=j0;
            b(i+1,2)=j1;
        end

        figure(1);
        plot(b(:,1),b(:,2));
        area_rand_train1(r,t)=0;
        for i=2:p0+q0+1,
            area_rand_train1(r,t)=area_rand_train1(r,t)+b(i,2)*(b(i,1)-b(i-1,1));
        end
        area_rand_train1(r,t);

% Total, average, geometric accuracy

        numtr0=[0 0];
        if sign0(k)==1
           for i=1:p0,
               if re0(i)<0
                  numtr0(1)=numtr0(1)+1;
               end
           end
           for i=p0+1:p0+q0,
               if re0(i)>0
                  numtr0(2)=numtr0(2)+1;
               end
           end
        end

        if sign0(k)==7
           numtr1=zeros(4,2);
           for i=1:p0,
               if (re(i,1)>0)&(re(i,2)<0)&(re(i,4)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(1,1)=numtr1(1,1)+1;
               end
               if (re(i,1)<0)&(re(i,3)>0)&(re(i,5)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(2,1)=numtr1(2,1)+1;
               end
               if (re(i,2)>0)&(re(i,6)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(3,1)=numtr1(3,1)+1;
               end
               if (re(i,3)<0)&(re(i,7)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(4,1)=numtr1(4,1)+1;
               end
           end

           for i=p0+1:p0+q0,
               if (re(i,1)>0)&(re(i,2)<0)&(re(i,4)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(1,2)=numtr1(1,2)+1;
               end
               if (re(i,1)<0)&(re(i,3)>0)&(re(i,5)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(2,2)=numtr1(2,2)+1;
               end
               if (re(i,2)>0)&(re(i,6)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(3,2)=numtr1(3,2)+1;
               end
               if (re(i,3)<0)&(re(i,7)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(4,2)=numtr1(4,2)+1;
               end
           end
        end

        if sign0(k)==12
           numtr1=zeros(8,2);
           for i=1:p0,
               if (re1(i,4)>0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,5)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(1,1)=numtr1(1,1)+1;
               end
               if (re1(i,4)>0)&((re1(i,1)<0)&(re1(i,3)>0)&(re1(i,6)<0))
                   numtr0(1)=numtr0(1)+1;
                   numtr1(2,1)=numtr1(2,1)+1;
               end
               if (re1(i,4)>0)&(re1(i,2)>0)&(re1(i,7)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(3,1)=numtr1(3,1)+1;
               end
               if (re1(i,4)>0)&(re1(i,3)<0)&(re1(i,8)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(4,1)=numtr1(4,1)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,9)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(5,1)=numtr1(5,1)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)<0)&(re1(i,3)>0)&(re1(i,10)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(6,1)=numtr1(6,1)+1;
               end
               if (re1(i,4)<0)&(re1(i,2)>0)&(re1(i,11)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(7,1)=numtr1(7,1)+1;
               end
               if (re1(i,4)<0)&(re1(i,3)<0)&(re1(i,12)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(8,1)=numtr1(8,1)+1;
               end
           end

           for i=p0+1:p0+q0,
               if (re1(i,4)>0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,5)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(1,2)=numtr1(1,2)+1;
               end
               if (re1(i,4)>0)&(re1(i,1)<0)&(re1(i,3)>0)&(re1(i,6)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(2,2)=numtr1(2,2)+1;
               end
               if (re1(i,4)>0)&(re1(i,2)>0)&(re1(i,7)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(3,2)=numtr1(3,2)+1;
               end
               if (re1(i,4)>0)&(re1(i,3)<0)&(re1(i,8)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(4,2)=numtr1(4,2)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,9)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(5,2)=numtr1(5,2)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)<0)&(re1(i,3)>0)&(re1(i,10)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(6,2)=numtr1(6,2)+1;
               end
               if (re1(i,4)<0)&(re1(i,2)>0)&(re1(i,11)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(7,2)=numtr1(7,2)+1;
               end
               if (re1(i,4)<0)&(re1(i,3)<0)&(re1(i,12)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(8,2)=numtr1(8,2)+1;
               end
           end
        end

        if sign0(k)==22
           numtr1=zeros(16,2);
           for i=1:p0,
               if re2(i,5)>0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,7)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(1,1)=numtr1(1,1)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,8)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(2,1)=numtr1(2,1)+1;
                  end
                  if re2(i,2)>0&re2(i,9)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(3,1)=numtr1(3,1)+1;
                  end
                  if re2(i,3)<0&re2(i,10)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(4,1)=numtr1(4,1)+1;
                  end
               end
               if re2(i,4)>0&re2(i,5)<0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,11)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(5,1)=numtr1(5,1)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,12)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(6,1)=numtr1(6,1)+1;
                  end
                  if re2(i,2)>0&re2(i,13)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(7,1)=numtr1(7,1)+1;
                  end
                  if re2(i,3)<0&re2(i,14)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(8,1)=numtr1(8,1)+1;
                  end
               end
               if re2(i,4)<0&re2(i,6)>0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,15)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(9,1)=numtr1(9,1)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,16)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(10,1)=numtr1(10,1)+1;
                  end
                  if re2(i,2)>0&re2(i,17)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(11,1)=numtr1(11,1)+1;
                  end
                  if re2(i,3)<0&re2(i,18)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(12,1)=numtr1(12,1)+1;
                  end
               end
               if re2(i,6)<0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,19)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(13,1)=numtr1(13,1)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,20)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(14,1)=numtr1(14,1)+1;
                  end
                  if re2(i,1)>0&re2(i,2)>0&re2(i,21)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(15,1)=numtr1(15,1)+1;
                  end
                  if re2(i,1)<0&re2(i,3)<0&re2(i,22)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(16,1)=numtr1(16,1)+1;
                  end
               end
           end

           for i=p0+1:p0+q0,
               if re2(i,5)>0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,7)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(1,2)=numtr1(1,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,8)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(2,2)=numtr1(2,2)+1;
                  end
                  if re2(i,1)>0&re2(i,2)>0&re2(i,9)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(3,2)=numtr1(3,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)<0&re2(i,10)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(4,2)=numtr1(4,2)+1;
                  end
               end
               if re2(i,4)>0&re2(i,5)<0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,11)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(5,2)=numtr1(5,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,12)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(6,2)=numtr1(6,2)+1;
                  end
                  if re2(i,1)>0&re2(i,2)>0&re2(i,13)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(7,2)=numtr1(7,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)<0&re2(i,14)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(8,2)=numtr1(8,2)+1;
                  end
               end
               if re2(i,4)<0&re2(i,6)>0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,15)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(9,2)=numtr1(9,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,16)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(10,2)=numtr1(10,2)+1;
                  end
                  if re2(i,1)>0&re2(i,2)>0&re2(i,17)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(11,2)=numtr1(11,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)<0&re2(i,18)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(12,2)=numtr1(12,2)+1;
                  end
               end
               if re2(i,6)<0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,19)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(13,2)=numtr1(13,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,20)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(14,2)=numtr1(14,2)+1;
                  end
                  if re2(i,1)>0&re2(i,2)>0&re2(i,21)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(15,2)=numtr1(15,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)<0&re2(i,22)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(16,2)=numtr1(16,2)+1;
                  end
               end
           end
        end
        numtr0
%        fprintf(fid,'%g  %g\n',numtr0);
%        if sign0(k)>1
%           numtr1
%           fprintf(fid,'%g  %g\n',numtr1);
%        end
%        fprintf(fid,'(r,t)=%d\n',r,t);

        err_rand_train1(r,t-1)=100*(1-(numtr0(1)+numtr0(2))/(p0+q0));
        err_rand_train1(r,c+t-2)=100*((p0-numtr0(1))/p0+(q0-numtr0(2))/q0)/2;
        err_rand_train1(r,2*c+t-3)=100*sqrt(((p0-numtr0(1))/p0)*((q0-numtr0(2))/q0));
        err_rand_train1(r,3*c+t-4)=100*area_rand_train1(r,t);
        err_rand_train1(r,4*c+t-5)=(err_rand_train1(r,t-1)+err_rand_train1(r,c+t-2)+err_rand_train1(r,2*c+t-3)+area_rand_train1(r,t)*100)/4;
    end
end

err_rand_train1
fid=fopen('D:\highimbalance\letreca\inv_dc1\result\err_rand_train1.txt','w');
    fprintf(fid,'%g %g %g %g %g\n',err_rand_train1');
fclose(fid);

t1=etime(clock,t0)

%clear all;
%t0=clock;
% Total, average, geometric accuracy; AUC:
% For the rand_test1 set.
load D:\highimbalance\letreca\data_5rand\rand_test1.txt;
load D:\highimbalance\letreca\nte.txt;
load D:\highimbalance\letreca\inv_dc1\result\w0.txt;
load D:\highimbalance\letreca\inv_dc1\result\w.txt;
load D:\highimbalance\letreca\inv_dc1\result\w1.txt;
load D:\highimbalance\letreca\inv_dc1\result\w2.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta0.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta1.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta2.txt;
load D:\highimbalance\letreca\inv_dc1\result\sign0.txt;

c=size(nte,1);
dim=size(rand_test1,2)-2;
%dim=8;
k=0;

for r=1:c-1,
    for t=r+1:c,
        k=k+1;
        p0=nte(r,3);
        q0=nte(t,3);
        a=rand_test1(nte(r,1):nte(r,2),1:dim);
        a(p0+1:p0+q0,:)=rand_test1(nte(t,1):nte(t,2),1:dim);
        re0=a*(w0(k,:))'+ones(p0+q0,1)*theta0(k);
        kk=sum(sign0(1:k));
        if sign0(k)==7
           re=a*(w(kk-6:kk,:))'+ones(p0+q0,1)*theta(kk-6:kk)';
        end
        if sign0(k)==12
           re1=a*(w1(kk-11:kk,:))'+ones(p0+q0,1)*theta1(kk-11:kk)';
        end
        if sign0(k)==22
           re2=a*(w2(kk-21:kk,:))'+ones(p0+q0,1)*theta2(kk-21:kk)';
        end

% AUC
        y=[];
        [y(:,1),y(:,2)]=sort(re0,'descend');

        b=[0 0];
        j10=0;
        j00=0;
        j1=0;
        j0=0;

        for i=1:p0+q0,
            if y(i,2)<=p0
               j10=j10+1;
               j1=j10/p0;
            else
               j00=j00+1;
               j0=j00/q0;
            end
            b(i+1,1)=j0;
            b(i+1,2)=j1;
        end

        figure(1);
        plot(b(:,1),b(:,2));
        area_rand_test1(r,t)=0;
        for i=2:p0+q0+1,
            area_rand_test1(r,t)=area_rand_test1(r,t)+b(i,2)*(b(i,1)-b(i-1,1));
        end
        area_rand_test1(r,t);

% Total, average, geometric accuracy

        numtr0=[0 0];
        if sign0(k)==1
           for i=1:p0,
               if re0(i)<0
                  numtr0(1)=numtr0(1)+1;
               end
           end
           for i=p0+1:p0+q0,
               if re0(i)>0
                  numtr0(2)=numtr0(2)+1;
               end
           end
        end

        if sign0(k)==7
           numtr1=zeros(4,2);
           for i=1:p0,
               if (re(i,1)>0)&(re(i,2)<0)&(re(i,4)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(1,1)=numtr1(1,1)+1;
               end
               if (re(i,1)<0)&(re(i,3)>0)&(re(i,5)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(2,1)=numtr1(2,1)+1;
               end
               if (re(i,2)>0)&(re(i,6)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(3,1)=numtr1(3,1)+1;
               end
               if (re(i,3)<0)&(re(i,7)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(4,1)=numtr1(4,1)+1;
               end
           end

           for i=p0+1:p0+q0,
               if (re(i,1)>0)&(re(i,2)<0)&(re(i,4)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(1,2)=numtr1(1,2)+1;
               end
               if (re(i,1)<0)&(re(i,3)>0)&(re(i,5)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(2,2)=numtr1(2,2)+1;
               end
               if (re(i,2)>0)&(re(i,6)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(3,2)=numtr1(3,2)+1;
               end
               if (re(i,3)<0)&(re(i,7)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(4,2)=numtr1(4,2)+1;
               end
           end
        end

        if sign0(k)==12
           numtr1=zeros(8,2);
           for i=1:p0,
               if (re1(i,4)>0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,5)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(1,1)=numtr1(1,1)+1;
               end
               if (re1(i,4)>0)&((re1(i,1)<0)&(re1(i,3)>0)&(re1(i,6)<0))
                   numtr0(1)=numtr0(1)+1;
                   numtr1(2,1)=numtr1(2,1)+1;
               end
               if (re1(i,4)>0)&(re1(i,2)>0)&(re1(i,7)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(3,1)=numtr1(3,1)+1;
               end
               if (re1(i,4)>0)&(re1(i,3)<0)&(re1(i,8)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(4,1)=numtr1(4,1)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,9)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(5,1)=numtr1(5,1)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)<0)&(re1(i,3)>0)&(re1(i,10)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(6,1)=numtr1(6,1)+1;
               end
               if (re1(i,4)<0)&(re1(i,2)>0)&(re1(i,11)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(7,1)=numtr1(7,1)+1;
               end
               if (re1(i,4)<0)&(re1(i,3)<0)&(re1(i,12)<0)
                   numtr0(1)=numtr0(1)+1;
                   numtr1(8,1)=numtr1(8,1)+1;
               end
           end

           for i=p0+1:p0+q0,
               if (re1(i,4)>0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,5)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(1,2)=numtr1(1,2)+1;
               end
               if (re1(i,4)>0)&(re1(i,1)<0)&(re1(i,3)>0)&(re1(i,6)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(2,2)=numtr1(2,2)+1;
               end
               if (re1(i,4)>0)&(re1(i,2)>0)&(re1(i,7)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(3,2)=numtr1(3,2)+1;
               end
               if (re1(i,4)>0)&(re1(i,3)<0)&(re1(i,8)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(4,2)=numtr1(4,2)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)>0)&(re1(i,2)<0)&(re1(i,9)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(5,2)=numtr1(5,2)+1;
               end
               if (re1(i,4)<0)&(re1(i,1)<0)&(re1(i,3)>0)&(re1(i,10)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(6,2)=numtr1(6,2)+1;
               end
               if (re1(i,4)<0)&(re1(i,2)>0)&(re1(i,11)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(7,2)=numtr1(7,2)+1;
               end
               if (re1(i,4)<0)&(re1(i,3)<0)&(re1(i,12)>0)
                   numtr0(2)=numtr0(2)+1;
                   numtr1(8,2)=numtr1(8,2)+1;
               end
           end
        end

        if sign0(k)==22
           numtr1=zeros(16,2);
           for i=1:p0,
               if re2(i,5)>0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,7)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(1,1)=numtr1(1,1)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,8)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(2,1)=numtr1(2,1)+1;
                  end
                  if re2(i,2)>0&re2(i,9)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(3,1)=numtr1(3,1)+1;
                  end
                  if re2(i,3)<0&re2(i,10)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(4,1)=numtr1(4,1)+1;
                  end
               end
               if re2(i,4)>0&re2(i,5)<0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,11)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(5,1)=numtr1(5,1)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,12)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(6,1)=numtr1(6,1)+1;
                  end
                  if re2(i,2)>0&re2(i,13)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(7,1)=numtr1(7,1)+1;
                  end
                  if re2(i,3)<0&re2(i,14)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(8,1)=numtr1(8,1)+1;
                  end
               end
               if re2(i,4)<0&re2(i,6)>0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,15)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(9,1)=numtr1(9,1)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,16)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(10,1)=numtr1(10,1)+1;
                  end
                  if re2(i,2)>0&re2(i,17)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(11,1)=numtr1(11,1)+1;
                  end
                  if re2(i,3)<0&re2(i,18)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(12,1)=numtr1(12,1)+1;
                  end
               end
               if re2(i,6)<0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,19)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(13,1)=numtr1(13,1)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,20)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(14,1)=numtr1(14,1)+1;
                  end
                  if re2(i,1)>0&re2(i,2)>0&re2(i,21)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(15,1)=numtr1(15,1)+1;
                  end
                  if re2(i,1)<0&re2(i,3)<0&re2(i,22)<0
                     numtr0(1)=numtr0(1)+1;
                     numtr1(16,1)=numtr1(16,1)+1;
                  end
               end
           end

           for i=p0+1:p0+q0,
               if re2(i,5)>0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,7)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(1,2)=numtr1(1,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,8)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(2,2)=numtr1(2,2)+1;
                  end
                  if re2(i,1)>0&re2(i,2)>0&re2(i,9)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(3,2)=numtr1(3,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)<0&re2(i,10)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(4,2)=numtr1(4,2)+1;
                  end
               end
               if re2(i,4)>0&re2(i,5)<0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,11)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(5,2)=numtr1(5,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,12)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(6,2)=numtr1(6,2)+1;
                  end
                  if re2(i,1)>0&re2(i,2)>0&re2(i,13)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(7,2)=numtr1(7,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)<0&re2(i,14)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(8,2)=numtr1(8,2)+1;
                  end
               end
               if re2(i,4)<0&re2(i,6)>0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,15)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(9,2)=numtr1(9,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,16)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(10,2)=numtr1(10,2)+1;
                  end
                  if re2(i,1)>0&re2(i,2)>0&re2(i,17)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(11,2)=numtr1(11,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)<0&re2(i,18)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(12,2)=numtr1(12,2)+1;
                  end
               end
               if re2(i,6)<0
                  if re2(i,1)>0&re2(i,2)<0&re2(i,19)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(13,2)=numtr1(13,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)>0&re2(i,20)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(14,2)=numtr1(14,2)+1;
                  end
                  if re2(i,1)>0&re2(i,2)>0&re2(i,21)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(15,2)=numtr1(15,2)+1;
                  end
                  if re2(i,1)<0&re2(i,3)<0&re2(i,22)>0
                     numtr0(2)=numtr0(2)+1;
                     numtr1(16,2)=numtr1(16,2)+1;
                  end
               end
           end
        end
        numtr0
%        fprintf(fid,'%g  %g\n',numtr0);
%        if sign0(k)>1
%           numtr1
%           fprintf(fid,'%g  %g\n',numtr1);
%        end
%        fprintf(fid,'(r,t)=%d\n',r,t);
        
        err_rand_test1(r,t-1)=100*(1-(numtr0(1)+numtr0(2))/(p0+q0));
        err_rand_test1(r,c+t-2)=100*((p0-numtr0(1))/p0+(q0-numtr0(2))/q0)/2;
        err_rand_test1(r,2*c+t-3)=100*sqrt(((p0-numtr0(1))/p0)*((q0-numtr0(2))/q0));
        err_rand_test1(r,3*c+t-4)=100*area_rand_test1(r,t);
        err_rand_test1(r,4*c+t-5)=(err_rand_test1(r,t-1)+err_rand_test1(r,c+t-2)+err_rand_test1(r,2*c+t-3)+area_rand_test1(r,t)*100)/4;
    end
end

err_rand_test1
fid=fopen('D:\highimbalance\letreca\inv_dc1\result\err_rand_test1.txt','w');
    fprintf(fid,'%g %g %g %g %g\n',err_rand_test1');
fclose(fid);

t1=etime(clock,t0)

%///////////////////////////////////////////////////////////////////////////%
%clear all;
%t0=clock;

% Majority vote for the rand_train1ing set

load D:\highimbalance\letreca\data_5rand\rand_train1.txt;
load D:\highimbalance\letreca\ntr.txt;
load D:\highimbalance\letreca\inv_dc1\result\w0.txt;
load D:\highimbalance\letreca\inv_dc1\result\w.txt;
load D:\highimbalance\letreca\inv_dc1\result\w1.txt;
load D:\highimbalance\letreca\inv_dc1\result\w2.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta0.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta1.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta2.txt;
load D:\highimbalance\letreca\inv_dc1\result\sign0.txt;

dim=size(w0,2);
nerr=[];
n0=size(rand_train1,1);
c=size(ntr,1);

n00=0;
for i=1:c-1,
    for j=i+1:c,
        n00=n00+1;
        sign00(i,j-1)=sign0(n00);
    end
end

sign000=sign00;
for i=1:c-1,
    for j=i+1:c,
        sign000(j,i)=-sign000(i,j-1);
    end
end

for i=1:n0,
    re0=rand_train1(i,1:dim)*w0'+theta0';
    re1=rand_train1(i,1:dim)*w'+theta';
    re2=rand_train1(i,1:dim)*w1'+theta1';
    re3=rand_train1(i,1:dim)*w2'+theta2';
    n00=0;
    for j=1:c-1,
        for k=j+1:c,
            n00=n00+1;
            a(j,k-1)=re0(n00);
        end
    end

    for j=1:c-1,
        for k=j+1:c,
            a(k,j)=-a(j,k-1);
        end
    end

    for q=1:c,
        k=c-1;
        for j=1:c-1,
            if abs(sign000(q,j))==1
               if a(q,j)<0
                  k=k-1;
               end
            end

            if sign000(q,j)==7
               re11=[];
               if q==1
                  kk=sum(sign00(1,1:j));
               else
                  kk=sum(sum(sign00(1:q-1,:)))+sum(sign00(q,1:j));
               end
               re11=re1(kk-6:kk);
               if ((re11(1)>0)&(re11(2)<0)&(re11(4)<0))|((re11(1)<0)&(re11(3)>0)&(re11(5)<0))|((re11(2)>0)&re11(6)<0)|((re11(3)<0)&re11(7)<0)
                  k=k-1;
               end
            end

            if sign000(q,j)==-7
               re11=[];
               j1=q-1;
               q1=j;
               if q1==1
                  kk=sum(sign00(1,1:j1));
               else
                  kk=sum(sum(sign00(1:q1-1,:)))+sum(sign00(q1,1:j1));
               end
               re11=re1(kk-6:kk);
               if ((re11(1)>0)&(re11(2)<0)&(re11(4)>0))|((re11(1)<0)&(re11(3)>0)&(re11(5)>0))|((re11(2)>0)&re11(6)>0)|((re11(3)<0)&re11(7)>0)
                  k=k-1;
               end
            end

            if sign000(q,j)==12
               re21=[];
               if q==1
                  kk=sum(sign00(1,1:j));
               else
                  kk=sum(sum(sign00(1:q-1,:)))+sum(sign00(q,1:j));
               end
               re21=re2(kk-11:kk);
               if ((re21(4)>0)&(((re21(1)>0)&(re21(2)<0)&(re21(5)<0))|((re21(1)<0)&(re21(3)>0)&(re21(6)<0))|((re21(2)>0)&(re21(7)<0))|((re21(3)<0)&(re21(8)<0))))|((re21(4)<0)&(((re21(1)>0)&(re21(2)<0)&(re21(9)<0))|((re21(1)<0)&(re21(3)>0)&(re21(10)<0))|((re21(2)>0)&(re21(11)<0))|((re21(3)<0)&(re21(12)<0))))
                  k=k-1;
               end
            end

            if sign000(q,j)==-12
               re21=[];
               j1=q-1;
               q1=j;
               if q1==1
                  kk=sum(sign00(1,1:j1));
               else
                  kk=sum(sum(sign00(1:q1-1,:)))+sum(sign00(q1,1:j1));
               end
               re21=re2(kk-11:kk);
               if ((re21(4)>0)&(((re21(1)>0)&(re21(2)<0)&(re21(5)>0))|((re21(1)<0)&(re21(3)>0)&(re21(6)>0))|((re21(2)>0)&(re21(7)>0))|((re21(3)<0)&(re21(8)>0))))|((re21(4)<0)&(((re21(1)>0)&(re21(2)<0)&(re21(9)>0))|((re21(1)<0)&(re21(3)>0)&(re21(10)>0))|((re21(2)>0)&(re21(11)>0))|((re21(3)<0)&(re21(12)>0))))
                  k=k-1;
               end
            end

            if sign000(q,j)==22
               re31=[];
               if q==1
                  kk=sum(sign00(1,1:j));
               else
                  kk=sum(sum(sign00(1:q-1,:)))+sum(sign00(q,1:j));
               end
               re31=re3(kk-21:kk);

               if re31(5)>0
                  if re31(1)>0&re31(2)<0&re31(7)<0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(8)<0
                     k=k-1;
                  end
                  if re31(2)>0&re31(9)<0
                     k=k-1;
                  end
                  if re31(3)<0&re31(10)<0
                     k=k-1;
                  end
               end
               if re31(4)>0&re31(5)<0
                  if re31(1)>0&re31(2)<0&re31(11)<0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(12)<0
                     k=k-1;
                  end
                  if re31(2)>0&re31(13)<0
                     k=k-1;
                  end
                  if re31(3)<0&re31(14)<0
                     k=k-1;
                  end
               end
               if re31(4)<0&re31(6)>0
                  if re31(1)>0&re31(2)<0&re31(15)<0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(16)<0
                     k=k-1;
                  end
                  if re31(2)>0&re31(17)<0
                     k=k-1;
                  end
                  if re31(3)<0&re31(18)<0
                     k=k-1;
                  end
               end
               if re31(6)<0
                  if re31(1)>0&re31(2)<0&re31(19)<0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(20)<0
                     k=k-1;
                  end
                  if re31(1)>0&re31(2)>0&re31(21)<0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)<0&re31(22)<0
                     k=k-1;
                  end
               end
            end

            if sign000(q,j)==-22
               re31=[];
               j1=q-1;
               q1=j;
               if q1==1
                  kk=sum(sign00(1,1:j1));
               else
                  kk=sum(sum(sign00(1:q1-1,:)))+sum(sign00(q1,1:j1));
               end
               re31=re3(kk-21:kk);

               if re31(5)>0
                  if re31(1)>0&re31(2)<0&re31(7)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(8)>0
                     k=k-1;
                  end
                  if re31(1)>0&re31(2)>0&re31(9)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)<0&re31(10)>0
                     k=k-1;
                  end
               end
               if re31(4)>0&re31(5)<0
                  if re31(1)>0&re31(2)<0&re31(11)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(12)>0
                     k=k-1;
                  end
                  if re31(1)>0&re31(2)>0&re31(13)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)<0&re31(14)>0
                     k=k-1;
                  end
               end
               if re31(4)<0&re31(6)>0
                  if re31(1)>0&re31(2)<0&re31(15)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(16)>0
                     k=k-1;
                  end
                  if re31(1)>0&re31(2)>0&re31(17)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)<0&re31(18)>0
                     k=k-1;
                  end
               end
               if re31(6)<0
                  if re31(1)>0&re31(2)<0&re31(19)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(20)>0
                     k=k-1;
                  end
                  if re31(1)>0&re31(2)>0&re31(21)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)<0&re31(22)>0
                     k=k-1;
                  end
               end
           end
        end
        result_tr(i,q)=k;
    end
end

nerr1=zeros(c,c);
nerr1_1=zeros(c,c);
kk0=0;
kk0_1=0;
for r=1:c,
    nerr(r)=0;
    nerr_1(r)=0;
    for i=ntr(r,1):ntr(r,2),
        b0=result_tr(i,r);
        b3=result_tr(i,:);
        b3(r)=-1;
        [amax,nmax]=max(b3);
        if amax>b0
           nerr(r)=nerr(r)+1;
           nerr1(r,nmax)=nerr1(r,nmax)+1;
           kk0=kk0+1;
        end
        if amax>=b0
           nerr_1(r)=nerr_1(r)+1;
           nerr1_1(r,nmax)=nerr1_1(r,nmax)+1;
           kk0_1=kk0_1+1;
        end
    end
end

nerr
nerr1
sum(nerr')
nerr_1
nerr1_1
sum(nerr_1')
1-sum(nerr')/ntr(c,2)
1-sum(nerr_1')/ntr(c,2)
fid=fopen('D:\highimbalance\letreca\inv_dc1\result\Number_error_rand_train1.txt','w');
    fprintf(fid,'\nFor the rand_train1 set\n');
    fprintf(fid,'\n%g %g\n',nerr');
    fprintf(fid,'\n%g %g\n',nerr1');
    fprintf(fid,'\n%g\n',sum(nerr'));
    fprintf(fid,'\n%g\n',1-sum(nerr')/ntr(c,2));
    fprintf(fid,'\n%g %g\n',nerr_1');
    fprintf(fid,'\n%g %g\n',nerr1_1');
    fprintf(fid,'\n%g\n',sum(nerr_1'));
    fprintf(fid,'\n%g\n',1-sum(nerr_1')/ntr(c,2));
fclose(fid);

t1=etime(clock,t0)

%///////////////////////////////////////////////////////////////////////////%
%clear all;
%t0=clock;
% Majority vote for the rand_test1 set

load D:\highimbalance\letreca\data_5rand\rand_test1.txt;
load D:\highimbalance\letreca\nte.txt;
load D:\highimbalance\letreca\inv_dc1\result\w0.txt;
load D:\highimbalance\letreca\inv_dc1\result\w.txt;
load D:\highimbalance\letreca\inv_dc1\result\w1.txt;
load D:\highimbalance\letreca\inv_dc1\result\w2.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta0.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta1.txt;
load D:\highimbalance\letreca\inv_dc1\result\theta2.txt;
load D:\highimbalance\letreca\inv_dc1\result\sign0.txt;

dim=size(w0,2);
nerr=[];
n0=size(rand_test1,1);
c=size(nte,1);

n00=0;
for i=1:c-1,
    for j=i+1:c,
        n00=n00+1;
        sign00(i,j-1)=sign0(n00);
    end
end

sign000=sign00;
for i=1:c-1,
    for j=i+1:c,
        sign000(j,i)=-sign000(i,j-1);
    end
end

for i=1:n0,
    re0=rand_test1(i,1:dim)*w0'+theta0';
    re1=rand_test1(i,1:dim)*w'+theta';
    re2=rand_test1(i,1:dim)*w1'+theta1';
    re3=rand_test1(i,1:dim)*w2'+theta2';
    n00=0;
    for j=1:c-1,
        for k=j+1:c,
            n00=n00+1;
            a(j,k-1)=re0(n00);
        end
    end

    for j=1:c-1,
        for k=j+1:c,
            a(k,j)=-a(j,k-1);
        end
    end

    for q=1:c,
        k=c-1;
        for j=1:c-1,
            if abs(sign000(q,j))==1
               if a(q,j)<0
                  k=k-1;
               end
            end

            if sign000(q,j)==7
               re11=[];
               if q==1
                  kk=sum(sign00(1,1:j));
               else
                  kk=sum(sum(sign00(1:q-1,:)))+sum(sign00(q,1:j));
               end
               re11=re1(kk-6:kk);
               if ((re11(1)>0)&(re11(2)<0)&(re11(4)<0))|((re11(1)<0)&(re11(3)>0)&(re11(5)<0))|((re11(2)>0)&re11(6)<0)|((re11(3)<0)&re11(7)<0)
                  k=k-1;
               end
            end

            if sign000(q,j)==-7
               re11=[];
               j1=q-1;
               q1=j;
               if q1==1
                  kk=sum(sign00(1,1:j1));
               else
                  kk=sum(sum(sign00(1:q1-1,:)))+sum(sign00(q1,1:j1));
               end
               re11=re1(kk-6:kk);
               if ((re11(1)>0)&(re11(2)<0)&(re11(4)>0))|((re11(1)<0)&(re11(3)>0)&(re11(5)>0))|((re11(2)>0)&re11(6)>0)|((re11(3)<0)&re11(7)>0)
                  k=k-1;
               end
            end

            if sign000(q,j)==12
               re21=[];
               if q==1
                  kk=sum(sign00(1,1:j));
               else
                  kk=sum(sum(sign00(1:q-1,:)))+sum(sign00(q,1:j));
               end
               re21=re2(kk-11:kk);
               if ((re21(4)>0)&(((re21(1)>0)&(re21(2)<0)&(re21(5)<0))|((re21(1)<0)&(re21(3)>0)&(re21(6)<0))|((re21(2)>0)&(re21(7)<0))|((re21(3)<0)&(re21(8)<0))))|((re21(4)<0)&(((re21(1)>0)&(re21(2)<0)&(re21(9)<0))|((re21(1)<0)&(re21(3)>0)&(re21(10)<0))|((re21(2)>0)&(re21(11)<0))|((re21(3)<0)&(re21(12)<0))))
                  k=k-1;
               end
            end

            if sign000(q,j)==-12
               re21=[];
               j1=q-1;
               q1=j;
               if q1==1
                  kk=sum(sign00(1,1:j1));
               else
                  kk=sum(sum(sign00(1:q1-1,:)))+sum(sign00(q1,1:j1));
               end
               re21=re2(kk-11:kk);
               if ((re21(4)>0)&(((re21(1)>0)&(re21(2)<0)&(re21(5)>0))|((re21(1)<0)&(re21(3)>0)&(re21(6)>0))|((re21(2)>0)&(re21(7)>0))|((re21(3)<0)&(re21(8)>0))))|((re21(4)<0)&(((re21(1)>0)&(re21(2)<0)&(re21(9)>0))|((re21(1)<0)&(re21(3)>0)&(re21(10)>0))|((re21(2)>0)&(re21(11)>0))|((re21(3)<0)&(re21(12)>0))))
                  k=k-1;
               end
            end

            if sign000(q,j)==22
               re31=[];
               if q==1
                  kk=sum(sign00(1,1:j));
               else
                  kk=sum(sum(sign00(1:q-1,:)))+sum(sign00(q,1:j));
               end
               re31=re3(kk-21:kk);

               if re31(5)>0
                  if re31(1)>0&re31(2)<0&re31(7)<0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(8)<0
                     k=k-1;
                  end
                  if re31(2)>0&re31(9)<0
                     k=k-1;
                  end
                  if re31(3)<0&re31(10)<0
                     k=k-1;
                  end
               end
               if re31(4)>0&re31(5)<0
                  if re31(1)>0&re31(2)<0&re31(11)<0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(12)<0
                     k=k-1;
                  end
                  if re31(2)>0&re31(13)<0
                     k=k-1;
                  end
                  if re31(3)<0&re31(14)<0
                     k=k-1;
                  end
               end
               if re31(4)<0&re31(6)>0
                  if re31(1)>0&re31(2)<0&re31(15)<0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(16)<0
                     k=k-1;
                  end
                  if re31(2)>0&re31(17)<0
                     k=k-1;
                  end
                  if re31(3)<0&re31(18)<0
                     k=k-1;
                  end
               end
               if re31(6)<0
                  if re31(1)>0&re31(2)<0&re31(19)<0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(20)<0
                     k=k-1;
                  end
                  if re31(1)>0&re31(2)>0&re31(21)<0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)<0&re31(22)<0
                     k=k-1;
                  end
               end
            end

            if sign000(q,j)==-22
               re31=[];
               j1=q-1;
               q1=j;
               if q1==1
                  kk=sum(sign00(1,1:j1));
               else
                  kk=sum(sum(sign00(1:q1-1,:)))+sum(sign00(q1,1:j1));
               end
               re31=re3(kk-21:kk);

               if re31(5)>0
                  if re31(1)>0&re31(2)<0&re31(7)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(8)>0
                     k=k-1;
                  end
                  if re31(1)>0&re31(2)>0&re31(9)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)<0&re31(10)>0
                     k=k-1;
                  end
               end
               if re31(4)>0&re31(5)<0
                  if re31(1)>0&re31(2)<0&re31(11)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(12)>0
                     k=k-1;
                  end
                  if re31(1)>0&re31(2)>0&re31(13)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)<0&re31(14)>0
                     k=k-1;
                  end
               end
               if re31(4)<0&re31(6)>0
                  if re31(1)>0&re31(2)<0&re31(15)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(16)>0
                     k=k-1;
                  end
                  if re31(1)>0&re31(2)>0&re31(17)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)<0&re31(18)>0
                     k=k-1;
                  end
               end
               if re31(6)<0
                  if re31(1)>0&re31(2)<0&re31(19)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)>0&re31(20)>0
                     k=k-1;
                  end
                  if re31(1)>0&re31(2)>0&re31(21)>0
                     k=k-1;
                  end
                  if re31(1)<0&re31(3)<0&re31(22)>0
                     k=k-1;
                  end
               end
           end
        end
        result_te(i,q)=k;
    end
end

nerr1=zeros(c,c);
nerr1_1=zeros(c,c);
kk0=0;
kk0_1=0;
for r=1:c,
    nerr(r)=0;
    nerr_1(r)=0;
    for i=nte(r,1):nte(r,2),
        b0=result_te(i,r);
        b3=result_te(i,:);
        b3(r)=-1;
        [amax,nmax]=max(b3);
        if amax>b0
           nerr(r)=nerr(r)+1;
           nerr1(r,nmax)=nerr1(r,nmax)+1;
           kk0=kk0+1;
        end
        if amax>=b0
           nerr_1(r)=nerr_1(r)+1;
           nerr1_1(r,nmax)=nerr1_1(r,nmax)+1;
           kk0_1=kk0_1+1;
        end
    end
end

nerr
nerr1
sum(nerr')
nerr_1
nerr1_1
sum(nerr_1')
1-sum(nerr')/nte(c,2)
1-sum(nerr_1')/nte(c,2)
fid=fopen('D:\highimbalance\letreca\inv_dc1\result\Number_error_rand_test1.txt','w');
    fprintf(fid,'\nFor the rand_test1 set\n');
    fprintf(fid,'\n%g %g\n',nerr');
    fprintf(fid,'\n%g %g\n',nerr1');
    fprintf(fid,'\n%g\n',sum(nerr'));
    fprintf(fid,'\n%g\n',1-sum(nerr')/nte(c,2));
    fprintf(fid,'\n%g %g\n',nerr_1');
    fprintf(fid,'\n%g %g\n',nerr1_1');
    fprintf(fid,'\n%g\n',sum(nerr_1'));
    fprintf(fid,'\n%g\n',1-sum(nerr_1')/nte(c,2));
fclose(fid);

t1=etime(clock,t0)