valuesOfX=linspace(0,pi/2,11);
valuesOfY=f(valuesOfX);
result1=trapezodialRule(valuesOfX,valuesOfY);
result2=simpsonRule(valuesOfX,valuesOfY);
actual=1;

actualError1= abs(1-result1);
actualError2= abs(1-result2);
fprintf("Result for trapezodial Rule:: %f\n",result1);
fprintf("Result for Simpson Rule: %f\n",result2);
fprintf("Actual error for trapezodial Rule: %f\n",actualError1);
fprintf("Actual error for Simpson Rule: %f\n",actualError2);

function result = trapezodialRule(valuesOfX,valuesOfY)
    sum=0;
    for i=2:11
        sum=sum+((valuesOfY(i)+valuesOfY(i-1))/2)*(valuesOfX(i)-valuesOfX(i-1));
    end
    result=sum;
end

function result=simpsonRule(valuesOfX,valuesOfY)
    sum=0;

    for i=2:2:size(valuesOfX,2)-1
        if(i~=2)
            sum=sum+4*valuesOfY(i)+2*valuesOfY(i-1);
        else
            sum=sum+4*valuesOfY(i);
        end
    end
    sum=sum+valuesOfY(1);
    sum=sum+valuesOfY(11);
    result = (pi/2)/(3*(size(valuesOfX,2)-1))*sum;
end

function value = f(x)
    value = sin(x);
end