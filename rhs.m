function [q] = rhs(A,b,u,N,jmax,stencil)

q = zeros(N,1);

switch stencil 
    case 9
        k=1;
        q(k) = b(k) - A(k,5)*u(k) - A(k,6)*u(k+1) - A(k,7)*u(k+jmax-1) - A(k,8)*u(k+jmax) - A(k,9)*u(k+jmax+1);
        for k=2:jmax-1
            q(k) = b(k) - A(k,4)*u(k-1) - A(k,5)*u(k) - A(k,6)*u(k+1) - A(k,7)*u(k+jmax-1) - A(k,8)*u(k+jmax) - A(k,9)*u(k+jmax+1);
        end
        k=jmax;
        q(k) = b(k) - A(k,3)*u(k-jmax+1) - A(k,4)*u(k-1) - A(k,5)*u(k) - A(k,6)*u(k+1) - A(k,7)*u(k+jmax-1) - A(k,8)*u(k+jmax) - A(k,9)*u(k+jmax+1);
        k=jmax+1;
        q(k) = b(k) - A(k,2)*u(k-jmax) - A(k,3)*u(k-jmax+1) - A(k,4)*u(k-1) - A(k,5)*u(k) - A(k,6)*u(k+1) - A(k,7)*u(k+jmax-1) - A(k,8)*u(k+jmax) - A(k,9)*u(k+jmax+1);
        for k = jmax+2:N-jmax-1
            q(k) = b(k) - A(k,1)*u(k-jmax-1) - A(k,2)*u(k-jmax) - A(k,3)*u(k-jmax+1) - A(k,4)*u(k-1) - A(k,5)*u(k) - A(k,6)*u(k+1) - A(k,7)*u(k+jmax-1) - A(k,8)*u(k+jmax) - A(k,9)*u(k+jmax+1);
        end
        k=N-jmax;
        q(k) = b(k) - A(k,1)*u(k-jmax-1) - A(k,2)*u(k-jmax) - A(k,3)*u(k-jmax+1) - A(k,4)*u(k-1) - A(k,5)*u(k) - A(k,6)*u(k+1) - A(k,7)*u(k+jmax-1) - A(k,8)*u(k+jmax);
        k = N-jmax+1;
        q(k) = b(k) - A(k,1)*u(k-jmax-1) - A(k,2)*u(k-jmax) - A(k,3)*u(k-jmax+1) - A(k,4)*u(k-1) - A(k,5)*u(k) - A(k,6)*u(k+1) - A(k,7)*u(k+jmax-1);
        for k=N-jmax+2:N-1
            q(k) = b(k) - A(k,1)*u(k-jmax-1) - A(k,2)*u(k-jmax) - A(k,3)*u(k-jmax+1) - A(k,4)*u(k-1) - A(k,5)*u(k) - A(k,6)*u(k+1);
        end
        k = N;
        q(k) = b(k) - A(k,1)*u(k-jmax-1) - A(k,2)*u(k-jmax) - A(k,3)*u(k-jmax+1) - A(k,4)*u(k-1) - A(k,5)*u(k);
        
    case 5
        k=1;
        q(k) = b(k) - A(k,3)*u(k) - A(k,4)*u(k+1) - A(k,5)*u(k+jmax);
        for k=2:jmax
            q(k) = b(k) - A(k,2)*u(k-1) - A(k,3)*u(k) - A(k,4)*u(k+1) - A(k,5)*u(k+jmax);
        end
        for k=jmax+1:N-jmax 
            q(k) = b(k) - A(k,1)*u(k-jmax) - A(k,2)*u(k-1) - A(k,3)*u(k) - A(k,4)*u(k+1) - A(k,5)*u(k+jmax); 
        end
        for k=N-jmax+1:N-1
            q(k) = b(k) - A(k,1)*u(k-jmax) - A(k,2)*u(k-1) - A(k,3)*u(k) - A(k,4)*u(k+1);
        end
        k=N;
        q(k) = b(k) - A(k,1)*u(k-jmax) - A(k,2)*u(k-1) - A(k,3)*u(k);    
end

end