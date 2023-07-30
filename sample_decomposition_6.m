% Sample decomposition 4

function  [a3,d3,n_tr3]=sample_decomposition_6(a0,w12,theta12,dim,p0,q0);
           result=[];
           result=a0(:,2:dim+1)*w12'+theta12;
           if p0==1
              mu=a0(p0,2:dim+1);
           else
              mu=mean(a0(1:p0,2:dim+1));
           end
           mu(2,:)=mean(a0(p0+1:p0+q0,2:dim+1));
           result_mu=mu*w12'+theta12;
           p11=0;
           p12=0;
           p13=0;
           p14=0;
           q11=0;
           q12=0;
           q13=0;
           q14=0;

           for i=1:p0,
               if result(i)>result_mu(1)
                   p13=p13+1;
                   a13(p13,:)=a0(i,:);
                   d13(p13,1)=1;
               elseif result(i)>=0
                   p11=p11+1;
                   a11(p11,:)=a0(i,:);
                   d11(p11,1)=1;
               elseif  result(i)>=result_mu(2)
                   p12=p12+1;
                   a12(p12,:)=a0(i,:);
                   d12(p12,1)=1;
               else
                   p14=p14+1;
                   a14(p14,:)=a0(i,:);
                   d14(p14,1)=1;
               end
           end

           for i=p0+1:p0+q0,
               if result(i)>result_mu(1)
                   q13=q13+1;
                   a13(p13+q13,:)=a0(i,:);
                   d13(p13+q13,1)=-1;
               elseif result(i)>=0
                   q11=q11+1;
                   a11(p11+q11,:)=a0(i,:);
                   d11(p11+q11,1)=-1;
               elseif  result(i)>=result_mu(2)
                   q12=q12+1;
                   a12(p12+q12,:)=a0(i,:);
                   d12(p12+q12,1)=-1;
               else
                   q14=q14+1;
                   a14(p14+q14,:)=a0(i,:);
                   d14(p14+q14,1)=-1;
               end
           end

           if p11>0|q11>0
              a3=a11;
              d3=d11;
           end
           if p12>0|q12>0
              a3(p11+q11+1:p11+q11+p12+q12,:)=a12;
              d3(p11+q11+1:p11+q11+p12+q12,:)=d12;
           end
           if p13>0|q13>0
              a3(p11+q11+p12+q12+1:p11+q11+p12+q12+p13+q13,:)=a13;
              d3(p11+q11+p12+q12+1:p11+q11+p12+q12+p13+q13,:)=d13;
           end
           if p14>0|q14>0
              a3(p11+q11+p12+q12+p13+q13+1:p11+q11+p12+q12+p13+q13+p14+q14,:)=a14;
              d3(p11+q11+p12+q12+p13+q13+1:p11+q11+p12+q12+p13+q13+p14+q14,:)=d14;
           end

           if p11>0
              n_tr3(1,1)=1;
              n_tr3(1,2)=p11;
           end
           n_tr3(1,3)=p11;
           n_tr3(1,4)=1;

           if q11>0
              n_tr3(2,1)=p11+1;
              n_tr3(2,2)=p11+q11;
           end
           n_tr3(2,3)=q11;
           n_tr3(2,4)=2;

           if p12>0
              n_tr3(3,1)=p11+q11+1;
              n_tr3(3,2)=p11+q11+p12;
           end
           n_tr3(3,3)=p12;
           n_tr3(3,4)=3;

           if q12>0
              n_tr3(4,1)=p11+q11+p12+1;
              n_tr3(4,2)=p11+q11+p12+q12;
           end
           n_tr3(4,3)=q12;
           n_tr3(4,4)=4;

           if p13>0
              n_tr3(5,1)=p11+q11+p12+q12+1;
              n_tr3(5,2)=p11+q11+p12+q12+p13;
           end
           n_tr3(5,3)=p13;
           n_tr3(5,4)=5;

           if q13>0
              n_tr3(6,1)=p11+q11+p12+q12+p13+1;
              n_tr3(6,2)=p11+q11+p12+q12+p13+q13;
           end
           n_tr3(6,3)=q13;
           n_tr3(6,4)=6;

           if p14>0
              n_tr3(7,1)=p11+q11+p12+q12+p13+q13+1;
              n_tr3(7,2)=p11+q11+p12+q12+p13+q13+p14;
           end
           n_tr3(7,3)=p14;
           n_tr3(7,4)=7;

           if q14>0
              n_tr3(8,1)=p11+q11+p12+q12+p13+q13+p14+1;
              n_tr3(8,2)=p11+q11+p12+q12+p13+q13+p14+q14;
           end
           n_tr3(8,3)=q14;
           n_tr3(8,4)=8;

%           fprintf('(p11, q11)+(p12, q12)+(p13, q13)+(p14, q14)=(%d, %d)+(%d, %d)+(%d, %d)+(%d, %d)\n',p11,q11,p12,q12,p13,q13,p14,q14);