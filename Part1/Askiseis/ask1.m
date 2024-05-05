

a1=-2;
b1=-1;

a2=-1;
b2=1;

a3=1.3; %sto 1 ksefevgei
b3=2;



printing(a1,b1);
printing(a2,b2);
printing(a3,b3);



function value = f(x)
    value = exp(sin(x).^3)+x.^6-2*x.^4-x.^3-1;
end

function value = fFirstDerivative(x)
    value = 3*exp(sin(x)^3)*cos(x)*sin(x)^2+6*x^5-8*x^3-3*x^2;
end

function [count,result] = bisec(a,b)
    err = (1/2)*10^(-5); 
    count=0;
    center = (a+b)/2;
    
    while(abs(b-a)>err)
        
        center = (a+b)/2;
        if (f(center)*f(a) < 0)
            b = center;
        else
            a=center;
        end
        count=count+1;
    end
    result=center;
       
    
end

function [count,result,quadratic] = newton(a,b)
    err = (1/2)*10^(-5);  
    quadratic=true;
    count=0;
    result=a+b/2;
    lastResult=9^100000;
    l=f(result)/fFirstDerivative(result);
    while (abs(result-lastResult) > err)
        lastResult=result;
        result=result-l;
        l=f(result)/fFirstDerivative(result);
        count=count+1;
    end
    if(l==0)
            quadratic=false;
    end
end


function [count,result] = secant(a,b)
    err = (1/2)*10^(-5); 
    count=0;
    
    while(abs(b-a)>err)
        count=count+1;
        temp=b-(f(b)*(b-a))/(f(b)-f(a));
        a=b;
        b=temp;

    end
    result=b;
end






function printing(a,b)
    fprintf("\n");
    fprintf("a=%.1f b=%.1f\n",a,b);
    fprintf("BISECTION:\n")
    if(f(a)*f(b)<0)
        [count,result]=bisec(a,b);
        fprintf("result = %f count = %d\n",result,count);
    else
        
        fprintf("Can not apply Bisection Method. (f(a)*f(b)>=0)\n") ;
    end

    fprintf("NEWTON-RAPHSON:\n")


    [count,result,quadratic]=newton(a,b);
    if(quadratic)
        fprintf("result = %f count = %d\n",result,count);
    else
        fprintf("f'(x)=0, not quadratic, nearest estimation was: %f, count=%d\n",result,count);
    end
    fprintf("SECANT METHOD:\n");
    [count,result]=secant(a,b);
    fprintf("result = %f count = %d\n",result,count);
end