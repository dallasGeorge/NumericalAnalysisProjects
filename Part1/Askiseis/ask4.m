n=15;
A=zeros(15);
q=0.15;

A(1,2)=1;
A(1,9)=1;

A(2,3)=1;
A(2,5)=1;
A(2,7)=1;

A(3,2)=1;
A(3,6)=1;
A(3,8)=1;

A(4,3)=1;
A(4,12)=1;

A(5,1)=1;
A(5,10)=1;

A(6,10)=1;
A(6,11)=1;

A(7,10)=1;
A(7,11)=1;

A(8,4)=1;
A(8,11)=1;

A(9,5)=1;
A(9,6)=1;
A(9,10)=1;

A(10,13)=1;

A(11,15)=1;

A(12,7)=1;
A(12,8)=1;
A(12,11)=1;

A(13,9)=1;
A(13,14)=1;

A(14,10)=1;
A(14,11)=1;
A(14,13)=1;
A(14,15)=1;

A(15,12)=1;
A(15,14)=1;

G=makeGoogle(A,n,q);

sumOfEachCol=zeros(1,n);
sumOfEachRow=zeros(1,n);
for j=1:n
    for i=1:n
        sumOfEachRow(1,j)=sumOfEachRow(1,j)+G(j,i);
        sumOfEachCol(1,j)=sumOfEachCol(1,j)+G(i,j);
    end
end

if(sumOfEachCol==ones(1,n))
    fprintf("matrix is left stochastic");
elseif(sumOfEachRow==ones(1,n))
    fprintf("matrix is right stochastic");
elseif(sumOfEachCol==ones(1,n) && sumOfEachRow==ones(1,n))
    fprintf("matrix is doubly stochastic");
else
    fprintf("not stochastic")
end

p=ones(n,1);
p = powerM(G,p);
fprintf("\np=\n");
for i=1:size(p,1)
    fprintf("p(%d,1) =%f\n",i,p(i,1));
end

%skopos veltiosh toy site 10
A2=A;
A2(13,14)=0;
A2(11,10)=1;
A2(8,10)=1;
A2(15,10)=1;
A2(13,10)=1;

G2=makeGoogle(A2,n,q);

p2=ones(n,1);
p2 = powerM(G2,p2);
fprintf("p2=\n");
for i=1:size(p2,1)
    fprintf("p2(%d,1) =%f\n",i,p2(i,1));
end

q2=0.02;
G3 = makeGoogle(A2,n,q2);
p3=ones(n,1);
p3 = powerM(G3,p3);
fprintf("p3=\n");
for i=1:size(p3,1)
    fprintf("p3(%d,1) =%f\n",i,p3(i,1));
end

q3=0.6;
G4 = makeGoogle(A2,n,q3);
p4=ones(n,1);
p4 = powerM(G4,p4);
fprintf("p4=\n");
for i=1:size(p4,1)
    fprintf("p4(%d,1) =%f\n",i,p4(i,1));
end

A3=A;
A3(8,11)=3;
A3(12,11)=3;
G5 = makeGoogle(A3,n,q);
p5=ones(n,1);
p5 = powerM(G5,p5);
fprintf("p5=\n");
for i=1:size(p5,1)
    fprintf("p5(%d,1) =%f\n",i,p5(i,1));
end

n2=14;
A4=A;
A4(10,:) = [];
A4(:,10) = [];
G6 = makeGoogle(A4,n2,q);
p6=ones(n2,1);
p6 = powerM(G6,p6);
fprintf("p6=\n");
for i=1:size(p6,1)
    fprintf("p6(%d,1) =%f\n",i,p6(i,1));
end


function p= powerM(G,p)
    eig1=p(1,1);
    diff=1;
    err=1/2*(10^-8);
    while diff>err
        p=G*p;
        eig2=p(1,1);
        p=p/eig2;
        diff=abs(eig1-eig2);
        eig1=eig2;
    end
    sum=0;
    for i=1:size(G,1)
        sum=sum+p(i,1);
    end
    p=p/sum;
end


function G= makeGoogle(A,n,q)
    G=zeros(n);
    for i=1:n
        for j=1:n
            nj=0;
            for j2=1:n
                nj=nj+A(j,j2);
            end
            G(i,j)=q/n+(A(j,i)*(1-q))/nj;
        end
    end
end
