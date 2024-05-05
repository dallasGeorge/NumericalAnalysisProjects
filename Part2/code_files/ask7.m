valuesOfX=(1:10);
valuesOfOPAP=[12.1300,12.1400,12.0100,12.2800,12.4800,12.7500,12.5900,12.1900,12.1500,12.0600];
valuesOfDogecoin=[0.060627,0.059288,0.060384,0.065962,0.064735,0.063444,0.062409,0.061678,0.062156,0.059513];
OPAPBirthday=12.1600;
OPAP5DaysAfterLastSample=12.7000;

DogecoinBirtday=0.060258;
Dogecoin5DaysAfterLastSample=0.058580;



result1a = leastSquare2(11,valuesOfX,valuesOfOPAP);
result1b = leastSquare2(11,valuesOfX,valuesOfDogecoin);
result2a = leastSquare3(11,valuesOfX,valuesOfOPAP);
result2b = leastSquare3(11,valuesOfX,valuesOfDogecoin);
result3a = leastSquare4(11,valuesOfX,valuesOfOPAP);
result3b = leastSquare4(11,valuesOfX,valuesOfDogecoin);

result4a = leastSquare2(15,valuesOfX,valuesOfOPAP);
result4b = leastSquare2(15,valuesOfX,valuesOfDogecoin);
result5a = leastSquare3(15,valuesOfX,valuesOfOPAP);
result5b = leastSquare3(15,valuesOfX,valuesOfDogecoin);
result6a = leastSquare4(15,valuesOfX,valuesOfOPAP);
result6b = leastSquare4(15,valuesOfX,valuesOfDogecoin);





fprintf("close values of OPAP for day of birthday: %f %f %f\n",result1a,result2a,result3a);
fprintf("close values of Dogecoin for day of birthday: %f %f %f\n",result1b,result2b,result3b);

fprintf("close values of OPAP for day of 5 days into future: %f %f %f\n",result4a,result5a,result6a);
fprintf("close values of Dogecoin for day of 5 days into future: %f %f %f\n",result4b,result5b,result6b);

function result = leastSquare2(x,valuesOfX,valuesOfY)
    A = zeros(10,3);
    
    for i=1:size(valuesOfX,2)

        A(i,3) =valuesOfX(i)^2;
        A(i,2) =valuesOfX(i);
        A(i,1) = 1;
        
    end
    ata=A'*A;

    b=A'*valuesOfY';
    abc=linsolve(ata,b);
    
    result=abc(1,1)+abc(2,1).*x+abc(3,1)*x.*x;

end

function result = leastSquare3(x,valuesOfX,valuesOfY)
    A = zeros(10,4);
    
    for i=1:size(valuesOfX,2)
        A(i,4)= valuesOfX(i)^3;
        A(i,3) =valuesOfX(i)^2;
        A(i,2) =valuesOfX(i);
        A(i,1) = 1;
        
    end
    ata=A'*A;

    b=A'*valuesOfY';
    abc=linsolve(ata,b);
    
    result=abc(1,1)+abc(2,1).*x+abc(3,1)*x.*x+abc(4,1).*x.*x.*x;

end

function result = leastSquare4(x,valuesOfX,valuesOfY)
    A = zeros(10,5);
    
    for i=1:size(valuesOfX,2)
        A(i,5)= valuesOfX(i)^4;
        A(i,4)= valuesOfX(i)^3;
        A(i,3) =valuesOfX(i)^2;
        A(i,2) =valuesOfX(i);
        A(i,1) = 1;
        
    end
    ata=A'*A;

    b=A'*valuesOfY';
    abc=linsolve(ata,b);
    
    result=abc(1,1)+abc(2,1).*x+abc(3,1)*x.*x+abc(4,1).*x.*x.*x+abc(5,1).*x.*x.*x.*x;

end