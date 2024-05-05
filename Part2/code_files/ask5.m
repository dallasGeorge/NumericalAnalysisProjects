valuesOfX=[-3.14, -2.53,-1.9, -1, 0, 0.51, 1.43, 1.84, 2.42,3.14];
x=linspace(-pi,pi,200);
result=LagrangeMethod(x,valuesOfX);
result2=LinearSplines(x,valuesOfX);
result3=leastSquare(x,valuesOfX);
plot(x,result);
grid on;
hold on;

plot(x,result2);
plot(x,result3);
for i=1:10
    plot(valuesOfX(i),f(valuesOfX(i)),'r*');
end
title("approximations of sin(x)","FontSize",8)
legend({'Lagrange','linear splines','least square'},'Location','southeast')
hold off;
figure(2) 
plot(x,abs(f(x)-result));
hold on;
plot(x,abs(f(x)-result2));
plot(x,abs(f(x)-result3));
title("Error of approximations compared to actual values of sin(x) for each x","FontSize",8)
legend({'Lagrange','linear splines','least square'},'Location','southeast')
hold off;
fprintf("Max error of lagrange: %f\n",max(abs(sin(x)-result)));
fprintf("Max error of linear splines: %f\n",max(abs(sin(x)-result2)));
fprintf("Max error of least square: %f\n",max(abs(sin(x)-result3)));

function value = f(x)
    value = round(sin(x),3);
end

function result = LagrangeMethod(x,valuesOfX)
    sum=0;
    for i=1:size(valuesOfX,2)
        a=1;
        b=1;
        for j=1:size(valuesOfX,2)
            if(i~=j)
                a=a.*(x-valuesOfX(j));
                b=b*(valuesOfX(i)-valuesOfX(j));
            end
        end
        sum=sum+f(valuesOfX(i))*a/b;
       
    end
    result=sum;
end

function result = LinearSplines(x,valuesOfX)
    
    
    for i=1:size(valuesOfX,2)-1
        result{i}= f(valuesOfX(i)) + (f(valuesOfX(i+1))-f(valuesOfX(i)))/(valuesOfX(i+1)-valuesOfX(i))*(x-valuesOfX(i));
    end
    k=1;
    for i=1:8
        point(k)=min(find(abs(result{i}-result{i+1})<=0.01));
        k=k+1;
    end
    point(9)=200;
    for i=2:8
        result{i}=result{i}(point(i-1):point(i)-1);
    end
    result{1}=result{1}(1:point(1)-1);
    result{9}=result{9}(point(8):point(9));
    result=cell2mat(result);

end

function result = leastSquare(x,valuesOfX)
    A = zeros(10,4);
    
    for i=1:size(valuesOfX,2)
     
        A(i,4)= valuesOfX(i)^3;
        A(i,3) =valuesOfX(i)^2;
        A(i,2) =valuesOfX(i);
        A(i,1) = 1;
        
    end
    ata=A'*A;

    b=A'*f(valuesOfX)';
    abc=linsolve(ata,b);
    
    result=abc(1,1)+abc(2,1).*x+abc(3,1).*x.^2+abc(4,1).*x.^3;

end