% Sample decomposition 3

function  [a51,d51,n_tr51]=sample_decomposition_7(a00,n_tr4,w5,theta5,dim);

           resultb=[];
           resultb=a00(:,2:dim+1)*w5'+theta5;
           p00=n_tr4(1,3);
           q00=n_tr4(2,3);

           n_tr51=zeros(4,4);
           a11=[];
           a21=[];
           a31=[];
           a41=[];
           d11=[];
           d21=[];
           d31=[];
           d41=[];           
           for i=1:p00,
               if resultb(i)>=0
                  n_tr51(1,3)=n_tr51(1,3)+1;
                  a11(n_tr51(1,3),:)=a00(i,:);
                  d11(n_tr51(1,3),1)=1;
               end

               if resultb(i)<0
                  n_tr51(3,3)=n_tr51(3,3)+1;
                  a31(n_tr51(3,3),:)=a00(i,:);
                  d31(n_tr51(3,3),1)=1;
               end
           end

           a21=[];
           a41=[];
           for i=p00+1:p00+q00,
               if resultb(i)>=0
                  n_tr51(2,3)=n_tr51(2,3)+1;
                  a21(n_tr51(2,3),:)=a00(i,:);
                  d21(n_tr51(2,3),1)=-1;
               end

               if resultb(i)<0
                  n_tr51(4,3)=n_tr51(4,3)+1;
                  a41(n_tr51(4,3),:)=a00(i,:);
                  d41(n_tr51(4,3),1)=-1;
               end
           end

           for i=1:4,
               n_tr51(i,4)=i;
           end
           
           if n_tr51(1,3)>0
              n_tr51(1,1)=1;
              n_tr51(1,2)=n_tr51(1,3);
              a51=a11;
              d51=d11;
           end

           if n_tr51(2,3)>0
              n_tr51(2,1)=n_tr51(1,3)+1;
              n_tr51(2,2)=sum(n_tr51(1:2,3));
              a51(n_tr51(1,3)+1:sum(n_tr51(1:2,3)),:)=a21;
              d51(n_tr51(1,3)+1:sum(n_tr51(1:2,3)),:)=d21;
           end

           if n_tr51(3,3)>0
              n_tr51(3,1)=sum(n_tr51(1:2,3))+1;
              n_tr51(3,2)=sum(n_tr51(1:3,3));
              a51(sum(n_tr51(1:2,3))+1:sum(n_tr51(1:3,3)),:)=a31;
              d51(sum(n_tr51(1:2,3))+1:sum(n_tr51(1:3,3)),:)=d31;
           end

%           if n_tr51(4,3)>0&(n_tr51(3,3)>0|n_tr51(2,3)>0|n_tr51(1,3)>0)
           if n_tr51(4,3)>0
              n_tr51(4,1)=sum(n_tr51(1:3,3))+1;
              n_tr51(4,2)=sum(n_tr51(1:4,3));
              a51(sum(n_tr51(1:3,3))+1:sum(n_tr51(1:4,3)),:)=a41;
              d51(sum(n_tr51(1:3,3))+1:sum(n_tr51(1:4,3)),:)=d41;
           end