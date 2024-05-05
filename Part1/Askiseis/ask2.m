a1=0;
b1=7/8;
extraForSecant1=(a1+b1)/2;

a2=0.9;
b2=2;
extraForSecant2=(a2+b2)/2;

a3=2.1;
b3=3;
extraForSecant3=(a3+b3)/2;


printing(a1,extraForSecant1,b1);
printing(a2,extraForSecant2,b2);
printing(a3,extraForSecant3,b3);

function value = f(x)
    value = 94*cos(x)^3-24*cos(x)+177*sin(x)^2-108*sin(x)^4-72*cos(x)^3*sin(x)^2-65;
end

function value = fFirstDerivative(x)
    value =(3.0/2.0)*(sin(x) + 6*sin(2*x)-15*sin(3*x))*(1-2*cos(x))^2;
end

function value = fSecondDerivative(x)
    value =(1-2*cos(x))*(6*sin(x)*(sin(x)-15*sin(3*x)+12*sin(x)*cos(x))-1.5*(2*cos(x)-1)*(cos(x)+12*cos(2*x)-45*cos(3*x)));
end
 
function [count,result] = bisec(a,b)
    err = (1/2.0)*10^(-5); 
    count=0;
    val=a;
    while(abs(b-a)>err)
        val = a + rand*(b-a);
        if (f(val)*f(a) < 0)
            b = val;
        else
            a=val;
        end
        count=count+1;
    end
    result=val;
       
    
end

function [count,result] = newton(a,b)
    err = (1/2.0)*10^(-5);  

    count=0;
    result=(a+b)/2;
    lastResult=a;
    l=1/(fFirstDerivative(result)/f(result)-(1/2.0)*(fSecondDerivative(result)/fFirstDerivative(result)));
    
    while (abs(result-lastResult) > err)
        lastResult=result;
        result=result-l;
        l=f(result)/fFirstDerivative(result);
        count=count+1;
    end

end


function [count,result] = secant(a,b,c)
    err = (1/2.0)*10^(-5); 
    count=0;
    while(abs(a-c)>err)
        q=f(a)/f(b);
        r=f(c)/f(b);
        s=f(c)/f(a);
        count=count+1;
        temp=c-((r*(r-q)*(c-b)+(1-r)*s*(c-a))/((q-1)*(r-1)*(s-1)));
        a=b;
        b=c;
        c=temp;

    end
    result=c;
end




function printing(a,extraForSecant,b)
    fprintf("\n");
    fprintf("a=%f b=%f extra point for secant=%f\n",a,b,extraForSecant);
    fprintf("BISECTION:\n")
    if(f(a)*f(b)<0)
        [count,result]=bisec(a,b);
        fprintf("result = %.7f count = %d\n",result,count);
    else
        
        fprintf("Can not apply Bisection Method. (f(a)*f(b)>=0)\n") ;
    end

    fprintf("NEWTON-RAPHSON:\n")
    [count,result]=newton(a,b);
    fprintf("result = %.7f count = %d\n",result,count);
    fprintf("SECANT METHOD:\n");
    [count,result]=secant(a,extraForSecant,b);
    fprintf("result = %.7f count = %d\n",result,count);
end