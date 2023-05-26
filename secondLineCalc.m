function [x_in2,y_in2,x_out2,y_out2] = secondLineCalc(imax,x_in,y_in,x_out,y_out,delta_in,delta_out,epsilon)

[x_in2,y_in2,x_out2,y_out2] = deal(zeros(imax,1));

%----------------------------second interal--------------------------------
i=1;
x_ksi = (x_in(2) - x_in(imax-1))/2;
y_ksi = (y_in(2) - y_in(imax-1))/2;
if (abs(y_ksi)<epsilon)
    x_in2(i) = x_in(i);
    if (y_in(i) > 0)
        y_in2(i) = y_in(i) + delta_in(i);
    else
        y_in2(i) = y_in(i) - delta_in(i);
    end
else
    x_in2(i) = x_in(i) + delta_in(i)/sqrt(1+(x_ksi/y_ksi)^2);
    y_in2(i) = y_in(i) - x_ksi/y_ksi*(x_in2(i)-x_in(i));
    x_eta = x_in2(i) - x_in(i);
    y_eta = y_in2(i) - y_in(i);
    J = x_ksi*y_eta - x_eta*y_ksi;
    if (J<=0)
        x_in2(i) = x_in(i) - delta_in(i)/sqrt(1+(x_ksi/y_ksi)^2);
        y_in2(i) = y_in(i) - x_ksi/y_ksi*(x_in2(i)-x_in(i));
    end
end

for i=2:imax-1
    x_ksi = (x_in(i+1) - x_in(i-1))/2;
    y_ksi = (y_in(i+1) - y_in(i-1))/2;
    if (abs(y_ksi)<epsilon)
        x_in2(i) = x_in(i);
        if (y_in(i) > 0)
            y_in2(i) = y_in(i) + delta_in(i);
        else
            y_in2(i) = y_in(i) - delta_in(i);
        end
    else
        x_in2(i) = x_in(i) + delta_in(i)/sqrt(1+(x_ksi/y_ksi)^2);
        y_in2(i) = y_in(i) - x_ksi/y_ksi*(x_in2(i)-x_in(i));
        x_eta = x_in2(i) - x_in(i);
        y_eta = y_in2(i) - y_in(i);
        J = x_ksi*y_eta - x_eta*y_ksi;
        if (J<=0)
            x_in2(i) = x_in(i) - delta_in(i)/sqrt(1+(x_ksi/y_ksi)^2);
            y_in2(i) = y_in(i) - x_ksi/y_ksi*(x_in2(i)-x_in(i));
        end
    end
end


i=imax;
x_ksi = (x_in(2) - x_in(imax-1))/2;
y_ksi = (y_in(2) - y_in(imax-1))/2;
if (abs(y_ksi)<epsilon)
    x_in2(i) = x_in(i);
    if (y_in(i) > 0)
        y_in2(i) = y_in(i) + delta_in(i);
    else
        y_in2(i) = y_in(i) - delta_in(i);
    end
else
    x_in2(i) = x_in(i) + delta_in(i)/sqrt(1+(x_ksi/y_ksi)^2);
    y_in2(i) = y_in(i) - x_ksi/y_ksi*(x_in2(i)-x_in(i));
    x_eta = x_in2(i) - x_in(i);
    y_eta = y_in2(i) - y_in(i);
    J = x_ksi*y_eta - x_eta*y_ksi;
    if (J<=0)
        x_in2(i) = x_in(i) - delta_in(i)/sqrt(1+(x_ksi/y_ksi)^2);
        y_in2(i) = y_in(i) - x_ksi/y_ksi*(x_in2(i)-x_in(i));
    end
end

%---------------------second external------------------------------------
i=1;
x_ksi = (x_out(2) - x_out(imax-1))/2;
y_ksi = (y_out(2) - y_out(imax-1))/2;
if (abs(y_ksi)<epsilon)
    x_out2(i) = x_out(i);
    if (y_out(i) > 0)
        y_out2(i) = y_out(i) - delta_out(i);
    else
        y_out2(i) = y_out(i) + delta_out(i);
    end
else
    x_out2(i) = x_out(i) + delta_out(i)/sqrt(1+(x_ksi/y_ksi)^2);
    y_out2(i) = y_out(i) + x_ksi/y_ksi*(x_out(i)-x_out2(i));
    x_eta = x_out(i) - x_out2(i);
    y_eta = y_out(i) - y_out2(i);
    J = x_ksi*y_eta - x_eta*y_ksi;
    if (J<=0)
        x_out2(i) = x_out(i) - delta_out(i)/sqrt(1+(x_ksi/y_ksi)^2);
        y_out2(i) = y_out(i) + x_ksi/y_ksi*(x_out(i)-x_out2(i));
    end
end

for i=2:imax-1
    x_ksi = (x_out(i+1) - x_out(i-1))/2;
    y_ksi = (y_out(i+1) - y_out(i-1))/2;
    if (abs(y_ksi)<epsilon)
        x_out2(i) = x_out(i);
        if (y_out(i) > 0)
            y_out2(i) = y_out(i) - delta_out(i);
        else
            y_out2(i) = y_out(i) + delta_out(i);
        end
    else
        x_out2(i) = x_out(i) + delta_out(i)/sqrt(1+(x_ksi/y_ksi)^2);
        y_out2(i) = y_out(i) + x_ksi/y_ksi*(x_out(i)-x_out2(i));
        x_eta = x_out(i) - x_out2(i);
        y_eta = y_out(i) - y_out2(i);
        J = x_ksi*y_eta - x_eta*y_ksi;
        if (J<=0)
            x_out2(i) = x_out(i) - delta_out(i)/sqrt(1+(x_ksi/y_ksi)^2);
            y_out2(i) = y_out(i) + x_ksi/y_ksi*(x_out(i)-x_out2(i));
        end
    end
end


i=imax;
x_ksi = (x_out(2) - x_out(imax-1))/2;
y_ksi = (y_out(2) - y_out(imax-1))/2;
if (abs(y_ksi)<epsilon)
    x_out2(i) = x_out(i);
    if (y_out(i) > 0)
        y_out2(i) = y_out(i) - delta_out(i);
    else
        y_out2(i) = y_out(i) + delta_out(i);
    end
else
    x_out2(i) = x_out(i) + delta_out(i)/sqrt(1+(x_ksi/y_ksi)^2);
    y_out2(i) = y_out(i) + x_ksi/y_ksi*(x_out(i)-x_out2(i));
    x_eta = x_out(i) - x_out2(i);
    y_eta = y_out(i) - y_out2(i);
    J = x_ksi*y_eta - x_eta*y_ksi;
    if (J<=0)
        x_out2(i) = x_out(i) - delta_out(i)/sqrt(1+(x_ksi/y_ksi)^2);
        y_out2(i) = y_out(i) + x_ksi/y_ksi*(x_out(i)-x_out2(i));
    end
end

end