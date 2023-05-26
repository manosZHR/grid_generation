%% Start

clear
close all
format long
clc

%% Input

imax=51;          %nodes in the ksi-direction
jmax=51;          %nodes in the eta-direction
maxiter_ex=1000;  %max iterations for the external loop
maxiter_in=500;  %max iterations for the internal loops
psi=0.6;            %damping coef. for the msip
tol_ex=1e-12;      %tolerance for the external loop
tol_in=1e-10;      %tolerance for the internal loops
omega_ex=0.6;       %relaxation factor for the external loop
omega_in=1;       %relaxation factor for the internal loop

r_in=1;           %radius of the inner boundary
r_out=10;         %radius of the outter boundary
delta_in=r_in/10*ones(imax,1);    %desired distance of the first internal grid line (internal boundary)
delta_out=r_out/5*ones(imax,1);   %desired distance of the first internal grid line (external boundary)
coef_in=1*ones(imax,1);           %coef. for the distance of the second internal grid line (internal boundary)
coef_out=0.5*ones(imax,1);          %coef. for the distance of the second internal grid line (external boundary)
epsilon=1e-14;    %tolerance for grid non-orthoginality
kappa_1=1;        %coef. for the interpolation of f1
kappa_2=10;        %coef. for the interpolation of f2


%% Data and matrices initialization

N=imax*jmax;      %total number of nodes

%boundary nodes and first internal grid line
[x_in,y_in,x_out,y_out]=cylinder(imax,r_in,r_out);
[f,x_in2,y_in2,x_out2,y_out2,x_in3,y_in3,x_out3,y_out3]=fCalc(x_in,y_in,x_out,y_out,imax,delta_in,delta_out,coef_in,coef_out,epsilon);

%% Init

%-----initialization-----%

[x,x_old,y,y_old]=deal(zeros(N,1));       %vectors of the new and old coordinates of each node
for i=1:imax
    for j=1:jmax
        k = (i-1)*jmax+j;
        x_old(k) = (r_in+(j-1)*(r_out-r_in)/(jmax-1))*cos((i-1)*(-2*pi/(imax-1)));
        y_old(k) = (r_in+(j-1)*(r_out-r_in)/(jmax-1))*sin((i-1)*(-2*pi/(imax-1)));
    end
end
x_old = 0.1*rand(N,1).*x_old;
y_old = 0.1*rand(N,1).*y_old;

figure
hold on
title('Initialization')
grid on
box on
axis tight
plot(x_old,y_old,'k.')

%% Solution

[dx,dy,w] = deal(zeros(N,1));

iter_ex=0;
figure
semilogy(0,1,'b.')
hold on
semilogy(0,1,'r.')
grid on
box on
axis tight
linkdata on
xlabel('$Iteration$','Interpreter','latex','FontSize',13)
ylabel('$Relative\,residual$','Interpreter','latex','FontSize',13)
while(iter_ex<maxiter_ex)
    iter_ex=iter_ex+1;
    [A,bx,by,metrics] = sysCon(x_old,y_old,imax,jmax,x_in,x_out,y_in,y_out,f,kappa_1,kappa_2,x_in2,y_in2,x_out2,y_out2,x_in3,y_in3,x_out3,y_out3);
    qx = rhs(A,bx,x_old,N,jmax,9);
    qy = rhs(A,by,y_old,N,jmax,9);
    sumx=0;
    sumy=0;
    for i=1:N
        sumx=sumx+qx(i)^2;
        sumy=sumy+qy(i)^2;
    end
    r0x=sqrt(sumx);
    r0y=sqrt(sumy);

    iter_in=0;
    x=x_old;
    y=y_old;
    LU = msip(A,psi,jmax,N,9);

    while (iter_in<maxiter_in)
        qx = rhs(A,bx,x,N,jmax,9);
        qy = rhs(A,by,y,N,jmax,9);
        sumx=0;
        sumy=0;
        for i=1:N
            sumx=sumx+qx(i)^2;
            sumy=sumy+qy(i)^2;
        end
        rx=sqrt(sumx)/r0x;
        ry=sqrt(sumy)/r0y;

        if ((rx<tol_in) && (ry<tol_in))
            break
        end
        iter_in=iter_in+1;

        %-----solution of the x-system-----%
        %forward
        k=1;
        w(k) = qx(k)/LU(k,5);
        for k=2:jmax-1
            w(k) = (qx(k) - LU(k,4)*w(k-1))/LU(k,5);
        end
        k=jmax;
        w(k) = (qx(k) - LU(k,3)*w(k-jmax+1) - LU(k,4)*w(k-1))/LU(k,5);
        k=jmax+1;
        w(k) = (qx(k) - LU(k,2)*w(k-jmax) - LU(k,3)*w(k-jmax+1) - LU(k,4)*w(k-1))/LU(k,5);
        for k=jmax+2:N
            w(k) = (qx(k) - LU(k,1)*w(k-jmax-1) - LU(k,2)*w(k-jmax) - LU(k,3)*w(k-jmax+1) - LU(k,4)*w(k-1))/LU(k,5);
        end
        %backward
        k=N;
        dx(k) = w(k);
        for k=N-1:-1:N-jmax+2
            dx(k) = w(k) - LU(k,6)*dx(k+1);
        end
        k=N-jmax+1;
        dx(k) = w(k) - LU(k,6)*dx(k+1) - LU(k,7)*dx(k+jmax-1);
        k = N-jmax;
        dx(k) = w(k) - LU(k,6)*dx(k+1) - LU(k,7)*dx(k+jmax-1) - LU(k,8)*dx(k+jmax);
        for k=N-jmax-1:-1:1
            dx(k) = w(k) - LU(k,6)*dx(k+1) - LU(k,7)*dx(k+jmax-1) - LU(k,8)*dx(k+jmax) - LU(k,9)*dx(k+jmax+1);
        end
        %x_new
        for i=1:N
            x(i) = omega_in*dx(i) + x(i);
        end

        %---solution of the y-system-----%
        %forward
        k=1;
        w(k) = qy(k)/LU(k,5);
        for k=2:jmax-1
            w(k) = (qy(k) - LU(k,4)*w(k-1))/LU(k,5);
        end
        k=jmax;
        w(k) = (qy(k) - LU(k,3)*w(k-jmax+1) - LU(k,4)*w(k-1))/LU(k,5);
        k=jmax+1;
        w(k) = (qy(k) - LU(k,2)*w(k-jmax) - LU(k,3)*w(k-jmax+1) - LU(k,4)*w(k-1))/LU(k,5);
        for k=jmax+2:N
            w(k) = (qy(k) - LU(k,1)*w(k-jmax-1) - LU(k,2)*w(k-jmax) - LU(k,3)*w(k-jmax+1) - LU(k,4)*w(k-1))/LU(k,5);
        end
        %backward
        k=N;
        dy(k) = w(k);
        for k=N-1:-1:N-jmax+2
            dy(k) = w(k) - LU(k,6)*dy(k+1);
        end
        k=N-jmax+1;
        dy(k) = w(k) - LU(k,6)*dy(k+1) - LU(k,7)*dy(k+jmax-1);
        k = N-jmax;
        dy(k) = w(k) - LU(k,6)*dy(k+1) - LU(k,7)*dy(k+jmax-1) - LU(k,8)*dy(k+jmax);
        for k=N-jmax-1:-1:1
            dy(k) = w(k) - LU(k,6)*dy(k+1) - LU(k,7)*dy(k+jmax-1) - LU(k,8)*dy(k+jmax) - LU(k,9)*dy(k+jmax+1);
        end
        %y_new
        for i=1:N
            y(i) = omega_in*dy(i) + y(i);
        end
    end

    sumx=0;
    sumy=0;
    sumx_old=0;
    sumy_old=0;
    for i=1:N
        sumx=sumx+(x(i)-x_old(i))^2;
        sumy=sumy+(y(i)-y_old(i))^2;
        sumx_old=sumx_old+x_old(i)^2;
        sumy_old=sumy_old+y_old(i)^2;
    end
    conv(iter_ex,1) = sqrt(sumx/sumx_old);
    conv(iter_ex,2) = sqrt(sumy/sumy_old);
    semilogy(iter_ex,conv(iter_ex,1),'b.')
    semilogy(iter_ex,conv(iter_ex,2),'r.')
    refresh
    drawnow
    fprintf('ex iter:\t%d\t\tin iter:\t%d\t\tresidual(x):\t%d\t\tresidual(y):\t%d\n',iter_ex,iter_in,conv(iter_ex,1),conv(iter_ex,2))

    if ((conv(iter_ex,1)<tol_ex) && (conv(iter_ex,2)<tol_ex))
        break
    else
        x_old = omega_ex*x + (1-omega_ex)*x_old;
        y_old = omega_ex*y + (1-omega_ex)*y_old;
    end
end   
legend('$x-residual$','$y-residual$','Interpreter','latex','FontSize',11)


%% Postprocess 

figure
hold on
grid on
box on
axis tight 
plot(x,y,'k.');

figure
hold on
box on
axis tight 
for j=1:jmax
    for i=1:imax-1
        k = (i-1)*jmax+j;
        plot([x(k) x(k+jmax)],[y(k) y(k+jmax)],'k','LineWidth',0.5)
    end
end
for i=1:imax
    for j=1:jmax-1
        k = (i-1)*jmax+j;
        plot([x(k) x(k+1)],[y(k) y(k+1)],'k','LineWidth',0.5)
    end
end
