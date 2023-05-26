function [f,x_in2,y_in2,x_out2,y_out2,x_in3,y_in3,x_out3,y_out3] = fCalc(x_in,y_in,x_out,y_out,imax,delta_in,delta_out,coef_in,coef_out,epsilon)

[x_in2,y_in2,x_out2,y_out2] = secondLineCalc(imax,x_in,y_in,x_out,y_out,delta_in,delta_out,epsilon);
[x_in3,y_in3,x_out3,y_out3] = secondLineCalc(imax,x_in2,y_in2,x_out2,y_out2,coef_in.*delta_in,coef_out.*delta_out,epsilon);

figure
hold on
grid on
box on
axis tight 
plot(x_in,y_in,'k.')
plot(x_in2,y_in2,'r.')
plot(x_in3,y_in3,'b.')
plot(x_out,y_out,'k.')
plot(x_out2,y_out2,'r.')
plot(x_out3,y_out3,'b.')

f = zeros(imax,4);

%-----f at eta=1 - internal boundary-----
i=1;
x_ksi = (x_in(2) - x_in(imax-1))/2;
y_ksi = (y_in(2) - y_in(imax-1))/2;
x_ksiksi = x_in(2) - 2*x_in(1) + x_in(imax-1);
y_ksiksi = y_in(2) - 2*y_in(1) + y_in(imax-1);
x_eta = x_in2(i) - x_in(i);
y_eta = y_in2(i) - y_in(i);
x_etaeta = (x_in3(i) - x_in(i))/2 - (x_in2(i) - x_in(i));
y_etaeta = (y_in3(i) - y_in(i))/2 - (y_in2(i) - y_in(i));
J = x_ksi*y_eta - x_eta*y_ksi;
invA = [y_eta -x_eta ; -y_ksi x_ksi]/J;
g_11 = x_ksi^2 + y_ksi^2;
g_22 = x_eta^2 + y_eta^2;
g_12 = x_ksi*x_eta + y_ksi*y_eta;
if (abs(g_12) >epsilon)
    fprintf("non-orthogonal grid at ksi = %d \n",i);
end
f(i,1:2) = (invA*[-g_22*x_ksiksi-g_11*x_etaeta; -g_22*y_ksiksi-g_11*y_etaeta]/J^2)';

for i=2:imax-1
    x_ksi = (x_in(i+1) - x_in(i-1))/2;
    y_ksi = (y_in(i+1) - y_in(i-1))/2;
    x_ksiksi = x_in(i+1) - 2*x_in(i) + x_in(i-1);
    y_ksiksi = y_in(i+1) - 2*y_in(i) + y_in(i-1);
    x_eta = x_in2(i) - x_in(i);
    y_eta = y_in2(i) - y_in(i);
    x_etaeta = (x_in3(i) - x_in(i))/2 - (x_in2(i) - x_in(i));
    y_etaeta = (y_in3(i) - y_in(i))/2 - (y_in2(i) - y_in(i));
    J = x_ksi*y_eta - x_eta*y_ksi;
    invA = [y_eta -x_eta ; -y_ksi x_ksi]/J;
    g_11 = x_ksi^2 + y_ksi^2;
    g_22 = x_eta^2 + y_eta^2;
    g_12 = x_ksi*x_eta + y_ksi*y_eta;
    if (abs(g_12) >epsilon)
        fprintf("non-orthogonal grid at ksi = %d \n",i);
    end
    f(i,1:2) = (invA*[-g_22*x_ksiksi-g_11*x_etaeta; -g_22*y_ksiksi-g_11*y_etaeta]/J^2)';
end

i = imax;
x_ksi = (x_in(2) - x_in(imax-1))/2;
y_ksi = (y_in(2) - y_in(imax-1))/2;
x_ksiksi = x_in(2) - 2*x_in(1) + x_in(imax-1);
y_ksiksi = y_in(2) - 2*y_in(1) + y_in(imax-1);
x_eta = x_in2(i) - x_in(i);
y_eta = y_in2(i) - y_in(i);
x_etaeta = (x_in3(i) - x_in(i))/2 - (x_in2(i) - x_in(i));
y_etaeta = (y_in3(i) - y_in(i))/2 - (y_in2(i) - y_in(i));
J = x_ksi*y_eta - x_eta*y_ksi;
invA = [y_eta -x_eta ; -y_ksi x_ksi]/J;
g_11 = x_ksi^2 + y_ksi^2;
g_22 = x_eta^2 + y_eta^2;
g_12 = x_ksi*x_eta + y_ksi*y_eta;
if (abs(g_12) >epsilon)
    fprintf("non-orthogonal grid at ksi = %d \n",i);
end
f(i,1:2) = (invA*[-g_22*x_ksiksi-g_11*x_etaeta; -g_22*y_ksiksi-g_11*y_etaeta]/J^2)';

%-----f at eta=etamax - external boundary-----
i = 1;
x_ksi = (x_out(2) - x_out(imax-1))/2;
y_ksi = (y_out(2) - y_out(imax-1))/2;
x_ksiksi = x_out(2) - 2*x_out(1) + x_out(imax-1);
y_ksiksi = y_out(2) - 2*y_out(1) + y_out(imax-1);
x_eta = x_out(i) - x_out2(i);
y_eta = y_out(i) - y_out2(i);
x_etaeta = (x_out(i) - x_out2(i)) - (x_out(i) - x_out3(i))/2;
J = x_ksi*y_eta - x_eta*y_ksi;
invA = [y_eta -x_eta ; -y_ksi x_ksi]/J;
g_11 = x_ksi^2 + y_ksi^2;
g_22 = x_eta^2 + y_eta^2;
g_12 = x_ksi*x_eta + y_ksi*y_eta;
if (abs(g_12) >epsilon)
    fprintf("non-orthogonal grid at ksi = %d \n",i);
end
f(i,3:4) = (invA*[-g_22*x_ksiksi-g_11*x_etaeta; -g_22*y_ksiksi-g_11*y_etaeta]/J^2)';
for i=2:imax-1
    x_ksi = (x_out(i+1) - x_out(i-1))/2;
    y_ksi = (y_out(i+1) - y_out(i-1))/2;
    x_ksiksi = x_out(i+1) - 2*x_out(i) + x_out(i-1);
    y_ksiksi = y_out(i+1) - 2*y_out(i) + y_out(i-1);
    x_eta = x_out(i) - x_out2(i);
    y_eta = y_out(i) - y_out2(i);
    x_etaeta = (x_out(i) - x_out2(i)) - (x_out(i) - x_out3(i))/2;
    J = x_ksi*y_eta - x_eta*y_ksi;
    invA = [y_eta -x_eta ; -y_ksi x_ksi]/J;
    g_11 = x_ksi^2 + y_ksi^2;
    g_22 = x_eta^2 + y_eta^2;
    g_12 = x_ksi*x_eta + y_ksi*y_eta;
    if (abs(g_12) >epsilon)
        fprintf("\n non-orthogonal grid at ksi = %d with g_12 = %f\n",i,g_12);
    end
    f(i,3:4) = (invA*[-g_22*x_ksiksi-g_11*x_etaeta; -g_22*y_ksiksi-g_11*y_etaeta]/J^2)';
end
i = imax;
x_ksi = (x_out(2) - x_out(imax-1))/2;
y_ksi = (y_out(2) - y_out(imax-1))/2;
x_ksiksi = x_out(2) - 2*x_out(1) + x_out(imax-1);
y_ksiksi = y_out(2) - 2*y_out(1) + y_out(imax-1);
x_eta = x_out(i) - x_out2(i);
y_eta = y_out(i) - y_out2(i);
x_etaeta = (x_out(i) - x_out2(i)) - (x_out(i) - x_out3(i))/2;
J = x_ksi*y_eta - x_eta*y_ksi;
invA = [y_eta -x_eta ; -y_ksi x_ksi]/J;
g_11 = x_ksi^2 + y_ksi^2;
g_22 = x_eta^2 + y_eta^2;
g_12 = x_ksi*x_eta + y_ksi*y_eta;
if (abs(g_12) >epsilon)
    fprintf("non-orthogonal grid at ksi = %d \n",i);
end
f(i,3:4) = (invA*[-g_22*x_ksiksi-g_11*x_etaeta; -g_22*y_ksiksi-g_11*y_etaeta]/J^2)';


end