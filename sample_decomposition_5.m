function  [a2,p1,q1,d1]=sample_decomposition_5(a1,w12,dim,d,p,q,p0,q0,fid);

           project00=[];
           project00=a1(:,2:dim+1)*w12';
           mu0=mean(project00(1:p));
           mu0(2)=mean(project00(p+1:p+q));
           a2=[];
           d1=[];

           if p/dim>8
              p1=0;
              for i=1:p,
                  if (project00(i)<mu0(1))&(project00(i)>mu0(2))
                      p1=p1+1;
                      a2(p1,:)=a1(i,:);
                      d1(p1,1)=d(i,1);
                  end
              end
           else
              p1=p;
              a2=a1(1:p1,:);
              d1=d(1:p1,1);
           end

           if q/dim>8
              q1=0;
              for i=p+1:p+q,
                  if (project00(i)<mu0(1))&(project00(i)>mu0(2))
                      q1=q1+1;
                      a2(p1+q1,:)=a1(i,:);
                      d1(p1+q1,1)=d(i,1);
                  end
              end
           else
              q1=q;
              a2(p1+1:p1+q1,:)=a1(p+1:p+q1,:);
              d1(p1+1:p1+q1,1)=d(p+1:p+q1,1);
           end

           fprintf('(p, q)=%d, %d\n',p1,q1);
           fprintf(fid,'(p, q)=%d, %d\n',p1,q1);