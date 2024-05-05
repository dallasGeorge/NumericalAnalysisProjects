
A=[5 2 2; 4 6 4; 4 2 2];
b=[4;-6;7];

MatrixForCholesky = [4 12 -16;12 37 -43; -16 -43 98];

x=luSolve(A,b);
fprintf("PA=LU result:")
display(x);

[A,b] = GaussSeidelMakeMatrix(10);
x=GaussSeidelSolve(A,b);
fprintf("gauss seidel result:");
display(x);
[A,b] = GaussSeidelMakeMatrix(10000);
x=GaussSeidelSolve(A,b);
fprintf("for n=10000 check the created x array in workspace of matlab.\n");

[L,check] = choleskyF(MatrixForCholesky);
if(check==1)
    fprintf("cholesky result");
    display(L);
else
    fprintf("Matrix is not Symmetric Positive Definite, cant apply Cholesky decomposotion.");
end    
function [x]=luSolve(A,b)
    [numOfRows,numOfCols] = size(A);
    
    P=eye(numOfRows); 
    L=P; 
    U=A; 

        %pivot
    for j=1:numOfCols
        absmax=0;
        maxrow=0;
        for i=j:numOfRows
            if(abs(U(i,j))>=absmax)
                absmax=abs(U(i,j));
                maxrow=i;
            end
        end
        U([j maxrow],:)=U([maxrow j],:);
        P([j maxrow],:)=P([maxrow j],:);
        for lin=1:j-1
            temp=L(j,lin);
            L(j,lin)=L(maxrow,lin);
            L(maxrow,lin)=temp;
                
        end
            
        piv=U(j,j);
  
        for n=j+1:numOfRows
            L(n,j)=U(n,j)/piv;
            for colu=j:numOfCols
                U(n,colu)=U(n,colu)-L(n,j)*U(j,colu);
            end
        end
            
    end

    
    y=fSub(L,P*b);
    x=bSub(U,y);


end
function x = fSub(L,b)
    
    m = size(L,1);
    x = zeros(m,1); %Matrix to be solved, y
    for j=1:m
        if (L(j,j) == 0)
         
            x(j,1) = 0;

            continue
        end

        temp = b(j,1);

        for i=1:j
            temp=temp - L(j,i) * x(i,1);
        end
    
        temp = temp / L(j,j);

        x(j,1) = temp;
    end
        

end
    

function x = bSub(U,b)
    m = size(U,1);
    x = zeros(m,1); 
    for j=m:-1:1 
        if (U(j,j) == 0)
    
            x(j,1) = 0;
            continue
        end
  
        temp = b(j,1);

        for i=j+1:m 
            temp = temp - U(j,i) * x(i,1);
        end
        temp = temp/ U(j,j);

        x(j,1) = temp;
   end
end



function [L,check]=choleskyF(A)
    if(issymmetric(A)==1 && all(eig(A)>0)==1)
        n=size(A,1);
        L=zeros(n);
        for k=1:n
            for i=1:k
                sum=0;
                sum2=0;
    
    
            
                if(k==i)
                    for j=1:k-1
                        sum=sum+L(k,j)^2;
                    end
                    L(k,i)=sqrt(A(k,i)-sum);
                else
                    for j=1:i-1
                        sum2=sum2+L(i,j)*L(k,j);
                    end
                    L(k,i)=(A(k,i)-sum2)/L(i,i);
                end
            end
        end
        check=1;
    else
        L=zeros(1);
        check=0;
    end
end


function [x]=GaussSeidelSolve(A,b)
    n=size(A,1);
    x=zeros(n,1);
    diff=99;
    lastx=x;
    while(diff>(1/2)*10^(-4))
        for i=1:n
            sum1=0;
            sum2=0;
            for j=1:i-1
                sum1=sum1+A(i,j)*x(j,1);
            end
            for j=i+1:n
                sum2=sum2+A(i,j)*lastx(j,1);
            end
            x(i,1)=1/A(i,i)*(b(i,1)-sum1-sum2);
        end
        diff=max(abs(x-lastx));
        lastx=x;
    end

end
function [A,b] = GaussSeidelMakeMatrix(n)
    A=zeros(n);
    b=ones(n,1);
    A(n,n)=5;
    for i=1:n-1
        A(i,i)=5;
        A(i+1,i)=-2;
        A(i,i+1)=-2;
    end
    b(1,1)=3;
    b(n,1)=3;
end

