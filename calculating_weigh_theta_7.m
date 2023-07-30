% calculating vertical weigh and theta %

function  [w5,theta5]=calculating_weigh_theta_7(mu01,mu02,w12,theta12,dim,fid);

           wa=w12;
           r_dim=dim;
           for i=1:dim,
               if abs(wa(i))<1.0e-15
                  r_dim=r_dim-1;
               end
           end

           for j=1:dim,
               wa(2,j)=j;
           end

           wa1=wa;
           for j=1:dim-1,
               for j1=j+1:dim,                   
                   if abs(wa1(1,j1))>abs(wa1(1,j))
                      wa0=wa1(:,j);
                      wa1(:,j)=wa1(:,j1);
                      wa1(:,j1)=wa0;
                   end
               end
           end

           wb(2,:)=wa1(2,:);
           wb(1,1)=1;
           for j=2:r_dim,
               rr=j-1;
               if abs(wa1(1,j))>=1.0e-15
                  wb(1,j)=((-1)^rr)*wa1(1,1)/wa1(1,j);
               end
           end

           if rem(r_dim,2)==1
              wb(1,r_dim)=0;
           end
          
           w51=wb;
           for j=1:dim-1,
               for j1=j+1:dim,
                   if w51(2,j1)<w51(2,j)
                      w50=w51(:,j);
                      w51(:,j)=w51(:,j1);
                      w51(:,j1)=w50;
                   end
               end
           end

           w5=w51(1,:);
           w5=w5/sqrt(w5*w5');
           if (mu01-mu02)*w5'<0
               w5=-w5;
           end
           t=-(w12*mu02'+theta12)/(w12*(mu01-mu02)');
           for i=1:dim,
               xx(i)=mu02(i)+(mu01(i)-mu02(i))*t;
           end
%xx
%//////////////////////////////////////////////////////////////////////////
%           for i=1:dim,
%               xx(i)=(mu02(i)*(w12*mu01')-mu01(i)*(w12*mu02')-theta12*(mu01(i)-mu02(i)))/(w12*(mu01-mu02)');
%           end
%/////////////////////////////////////////////////////////////////////////
%xx
           theta5=-w5*xx';

%           w12
%           w5
%           w12*w5'
           fprintf('\n w12*w5=%g\n',w12*w5');
           fprintf(fid,'\n w12*w5=%g\n',w12*w5');