% Sample decomposition 4

function  [a71,d71,n_tr71]=sample_decomposition_8(a00,n_tr6,w5,theta5,xigma3,dim);

           resultb=[];
           resultb=a00(:,2:dim+1)*w5'+theta5;
           p0=n_tr6(1,3);
           q0=n_tr6(2,3);

           n_tr71=zeros(4,4);
           a11=[];
           a31=[];
           for i=1:p0,
               if resultb(i)>=xigma3
                  n_tr71(1,3)=n_tr71(1,3)+1;
                  a11(n_tr71(1,3),:)=a00(i,:);
                  d11(n_tr71(1,3),1)=1;
               end

               if resultb(i)<xigma3
                  n_tr71(3,3)=n_tr71(3,3)+1;
                  a31(n_tr71(3,3),:)=a00(i,:);
                  d31(n_tr71(3,3),1)=1;
               end
           end

           a21=[];
           a41=[];
           for i=p0+1:p0+q0,
               if resultb(i)>=xigma3
                  n_tr71(2,3)=n_tr71(2,3)+1;
                  a21(n_tr71(2,3),:)=a00(i,:);
                  d21(n_tr71(2,3),1)=-1;
               end

               if resultb(i)<xigma3
                  n_tr71(4,3)=n_tr71(4,3)+1;
                  a41(n_tr71(4,3),:)=a00(i,:);
                  d41(n_tr71(4,3),1)=-1;
               end
           end

           for i=1:4,
               n_tr71(i,4)=i;
           end
           
           if n_tr71(1,3)>0
              n_tr71(1,1)=1;
              n_tr71(1,2)=n_tr71(1,3);
              a71=a11;
              d71=d11;
           end
           if n_tr71(2,3)>0
              n_tr71(2,1)=n_tr71(1,3)+1;
              n_tr71(2,2)=n_tr71(1,3)+n_tr71(2,3);
              a71(n_tr71(1,3)+1:n_tr71(1,3)+n_tr71(2,3),:)=a21;
              d71(n_tr71(1,3)+1:n_tr71(1,3)+n_tr71(2,3),:)=d21;
           end
           if n_tr71(3,3)>0
              n_tr71(3,1)=n_tr71(1,3)+n_tr71(2,3)+1;
              n_tr71(3,2)=n_tr71(1,3)+n_tr71(2,3)+n_tr71(3,3);
              a71(n_tr71(1,3)+n_tr71(2,3)+1:n_tr71(1,3)+n_tr71(2,3)+n_tr71(3,3),:)=a31;
              d71(n_tr71(1,3)+n_tr71(2,3)+1:n_tr71(1,3)+n_tr71(2,3)+n_tr71(3,3),:)=d31;
           end
           if n_tr71(4,3)>0
              n_tr71(4,1)=n_tr71(1,3)+n_tr71(2,3)+n_tr71(3,3)+1;
              n_tr71(4,2)=n_tr71(1,3)+n_tr71(2,3)+n_tr71(3,3)+n_tr71(4,3);
              a71(n_tr71(1,3)+n_tr71(2,3)+n_tr71(3,3)+1:n_tr71(1,3)+n_tr71(2,3)+n_tr71(3,3)+n_tr71(4,3),:)=a41;
              d71(n_tr71(1,3)+n_tr71(2,3)+n_tr71(3,3)+1:n_tr71(1,3)+n_tr71(2,3)+n_tr71(3,3)+n_tr71(4,3),:)=d41;
           end
%           n_tr71


