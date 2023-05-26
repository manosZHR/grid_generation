function [A,bx,by,metrics] = sysCon(x,y,imax,jmax,x_in,x_out,y_in,y_out,f,kappa_1,kappa_2,x_in2,y_in2,x_out2,y_out2,x_in3,y_in3,x_out3,y_out3)

N = imax*jmax;
metrics=zeros(N,5);   %metrics of each node
A=zeros(N,9);         %coefficient matrix
bx=zeros(N,1);        %vector - x
by=zeros(N,1);        %vector - y

%------------metrics------------%
for i=2:imax-1
    for j=2:jmax-1
        k=(i-1)*jmax+j;
        metrics(k,1) = (x(k+jmax)-x(k-jmax))/2;                                 %x_ksi
        metrics(k,2) = (x(k+1)-x(k-1))/2;                                       %x_eta
        metrics(k,3) = (y(k+jmax)-y(k-jmax))/2;                                 %y_ksi
        metrics(k,4) = (y(k+1)-y(k-1))/2;                                       %y_eta
        metrics(k,5) = metrics(k,1)*metrics(k,4) - metrics(k,2)*metrics(k,3);   %J
    end
end
%split line (1)
i=1;
j=1;
k=(i-1)*jmax+j;
metrics(k,1) = x(k+jmax)-x(k);
metrics(k,2) = x(k+1)-x(k);
metrics(k,3) = y(k+jmax)-y(k);
metrics(k,4) = y(k+1)-y(k);
metrics(k,5) = metrics(k,1)*metrics(k,4) - metrics(k,2)*metrics(k,3);
for j=2:jmax-1
    k=(i-1)*jmax+j;
    metrics(k,1) = x(k+jmax)-x(k);
    metrics(k,2) = (x(k+1)-x(k-1))/2;
    metrics(k,3) = y(k+jmax)-y(k);
    metrics(k,4) = (y(k+1)-y(k-1))/2;
    metrics(k,5) = metrics(k,1)*metrics(k,4) - metrics(k,2)*metrics(k,3);
end
j=jmax;
k=(i-1)*jmax+j;
metrics(k,1) = x(k+jmax)-x(k);
metrics(k,2) = x(k)-x(k-1);
metrics(k,3) = y(k+jmax)-y(k);
metrics(k,4) = y(k)-y(k-1);
metrics(k,5) = metrics(k,1)*metrics(k,4) - metrics(k,2)*metrics(k,3);
%split line (max)
i=imax;
j=1;
k=(i-1)*jmax+j;
metrics(k,1) = x(k)-x(k-jmax);
metrics(k,2) = x(k+1)-x(k);
metrics(k,3) = y(k)-y(k-jmax);
metrics(k,4) = y(k+1)-y(k);
metrics(k,5) = metrics(k,1)*metrics(k,4) - metrics(k,2)*metrics(k,3);
for j=2:jmax-1
    k=(i-1)*jmax+j;
    metrics(k,1) = x(k)-x(k-jmax);
    metrics(k,2) = (x(k+1)-x(k-1))/2;
    metrics(k,3) = y(k)-y(k-jmax);
    metrics(k,4) = (y(k+1)-y(k-1))/2;
    metrics(k,5) = metrics(k,1)*metrics(k,4) - metrics(k,2)*metrics(k,3);
end
j=jmax;
k=(i-1)*jmax+j;
metrics(k,1) = x(k)-x(k-jmax);
metrics(k,2) = x(k)-x(k-1);
metrics(k,3) = y(k)-y(k-jmax);
metrics(k,4) = y(k)-y(k-1);
metrics(k,5) = metrics(k,1)*metrics(k,4) - metrics(k,2)*metrics(k,3);
%inner boundary
j=1;
for i=2:imax-1
    k=(i-1)*jmax+j;
    metrics(k,1) = (x(k+jmax)-x(k-jmax))/2;
    metrics(k,2) = x(k+1)-x(k);
    metrics(k,3) = (y(k+jmax)-y(k-jmax))/2;
    metrics(k,4) = y(k+1)-y(k);
    metrics(k,5) = metrics(k,1)*metrics(k,4) - metrics(k,2)*metrics(k,3);
end
%outter boundary
j=jmax;
for i=2:imax-1
    k=(i-1)*jmax+j;
    metrics(k,1) = (x(k+jmax)-x(k-jmax))/2;
    metrics(k,2) = x(k)-x(k-1);
    metrics(k,3) = (y(k+jmax)-y(k-jmax))/2;
    metrics(k,4) = y(k)-y(k-1);
    metrics(k,5) = metrics(k,1)*metrics(k,4) - metrics(k,2)*metrics(k,3);
end

%-----matrices construction-----%
%inner boundary
j=1;
for i=1:imax
    k = (i-1)*jmax+j;
    A(k,5) = 1;
    bx(k) = x_in(i);
    by(k) = y_in(i);
end
%outter boundary
j=jmax;
for i=1:imax
    k = (i-1)*jmax+j;
    A(k,5) = 1;
    bx(k) = x_out(i);
    by(k) = y_out(i);
end
%split line (1)
i=1;
for j=2:jmax-1
    k=(i-1)*jmax+j;
    g_11=metrics(k,1)^2+metrics(k,3)^2;
    g_12=metrics(k,1)*metrics(k,2) + metrics(k,3)*metrics(k,4);
    g_22=metrics(k,2)^2+metrics(k,4)^2;
    J=metrics(k,5);
    f1=f(i,1)*exp(-kappa_1*(j-1))+f(i,3)*exp(-kappa_1*(jmax-j));
    f2=f(i,2)*exp(-kappa_2*(j-1))+f(i,4)*exp(-kappa_2*(jmax-j));
    A(k,4)=g_11-J^2*f2/2;
    A(k,5)=-2*g_11-2*g_22;
    A(k,6)=g_11+J^2*f2/2;
    A(k,7)=g_12/2;
    A(k,8)=g_22+J^2*f1/2;
    A(k,9)=-g_12/2;
    bx(k)=-(-g_12/2)*x((imax-2)*jmax+j-1)-(g_22-J^2*f1/2)*x((imax-2)*jmax+j)-(g_12/2)*x((imax-2)*jmax+j+1);
    by(k)=-(-g_12/2)*y((imax-2)*jmax+j-1)-(g_22-J^2*f1/2)*y((imax-2)*jmax+j)-(g_12/2)*y((imax-2)*jmax+j+1);
end
%split line (max)
i=imax;
for j=2:jmax-1
    k=(i-1)*jmax+j;
    g_11=metrics(k,1)^2+metrics(k,3)^2;
    g_12=metrics(k,1)*metrics(k,2) + metrics(k,3)*metrics(k,4);
    g_22=metrics(k,2)^2+metrics(k,4)^2;
    J=metrics(k,5);
    f1=f(i,1)*exp(-kappa_1*(j-1))+f(i,3)*exp(-kappa_1*(jmax-j));
    f2=f(i,2)*exp(-kappa_2*(j-1))+f(i,4)*exp(-kappa_2*(jmax-j));
    A(k,1)=-g_12/2;
    A(k,2)=g_22-J^2*f1/2;
    A(k,3)=g_12/2;
    A(k,4)=g_11-J^2*f2/2;
    A(k,5)=-2*g_11-2*g_22;
    A(k,6)=g_11+J^2*f2/2;
    bx(k)=-(-g_12/2)*x(jmax+j+1)-(g_22+J^2*f1/2)*x(jmax+j)-(g_12/2)*x(jmax+j-1);
    by(k)=-(-g_12/2)*y(jmax+j+1)-(g_22+J^2*f1/2)*y(jmax+j)-(g_12/2)*y(jmax+j-1);
end
%inner nodes
for i=2:imax-1
    for j=2:jmax-1
        k=(i-1)*jmax+j;
        g_11=metrics(k,1)^2+metrics(k,3)^2;
        g_12=metrics(k,1)*metrics(k,2) + metrics(k,3)*metrics(k,4);
        g_22=metrics(k,2)^2+metrics(k,4)^2;
        J=metrics(k,5);
        f1=f(i,1)*exp(-kappa_1*(j-1))+f(i,3)*exp(-kappa_1*(jmax-j));
        f2=f(i,2)*exp(-kappa_2*(j-1))+f(i,4)*exp(-kappa_2*(jmax-j));
        A(k,1)=-g_12/2;
        A(k,2)=g_22-J^2*f1/2;
        A(k,3)=g_12/2;
        A(k,4)=g_11-J^2*f2/2;
        A(k,5)=-2*g_11-2*g_22;
        A(k,6)=g_11+J^2*f2/2;
        A(k,7)=g_12/2;
        A(k,8)=g_22+J^2*f1/2;
        A(k,9)=-g_12/2;
        bx(k)=0;
        by(k)=0;
    end
end

% extra boundary conditions
% j=2;
% for i=1:imax
%     k = (i-1)*jmax+j;
%     A(k,:) = zeros(1,9);
%     A(k,5) = 1;
%     bx(k) = x_in2(i);
%     by(k) = y_in2(i);
% end
% j=3;
% for i=1:imax
%     k = (i-1)*jmax+j;
%     A(k,:) = zeros(1,9);
%     A(k,5) = 1;
%     bx(k) = x_in3(i);
%     by(k) = y_in3(i);
% end
% %outter boundary
% j=jmax-1;
% for i=1:imax
%     k = (i-1)*jmax+j;
%     A(k,:) = zeros(1,9);
%     A(k,5) = 1;
%     bx(k) = x_out2(i);
%     by(k) = y_out2(i);
% end
% j=jmax-2;
% for i=1:imax
%     k = (i-1)*jmax+j;
%     A(k,:) = zeros(1,9);
%     A(k,5) = 1;
%     bx(k) = x_out3(i);
%     by(k) = y_out3(i);
% end

end