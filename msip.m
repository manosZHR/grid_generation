function LU = msip(A,psi,jmax,N,stencil)

switch stencil
    
    case 5
        
        LU = zeros(N,stencil);
        
        k=1;
        %A(k,3) = A(k,3)            %e(k) = E(k)
        A(k,4) = A(k,4)/A(k,3);     %f(k) = F(k)/e(k)
        A(k,5) = A(k,5)/A(k,3);     %h(k) = H(k)/e(k)
        
        for k=2:jmax
            A(k,2) = A(k,2)/(1+psi*A(k-1,5));   %d(k) = D(k)/[1+psi*h(k-1)]
            A(k,3) = A(k,3) - A(k,2)*A(k-1,4) + psi*A(k,2)*A(k-1,5);
            A(k,4) = A(k,4)/A(k,3);     %f(k) = F(k)/e(k)
            A(k,5) = (A(k,5) - psi*A(k,2)*A(k-1,5))/A(k,3);
        end
        
        for k=jmax+1:N-jmax
            A(k,1) = A(k,1)/(1+psi*A(k-jmax,4));
            A(k,2) = A(k,2)/(1+psi*A(k-1,5)); 
            A(k,3) = A(k,3) - A(k,2)*A(k-1,4) - A(k,1)*A(k-jmax,5) + psi*(A(k,2)*A(k-1,5)+ A(k,1)*A(k-jmax,4));
            A(k,4) = (A(k,4) - psi*A(k,1)*A(k-jmax,4))/A(k,3);
            A(k,5) = (A(k,5) - psi*A(k,2)*A(k-1,5))/A(k,3);
        end
        
        for k=N-jmax+1:N-1
            A(k,1) = A(k,1)/(1+psi*A(k-jmax,4));
            A(k,2) = A(k,2)/(1+psi*A(k-1,5));
            A(k,3) = A(k,3) - A(k,2)*A(k-1,4) - A(k,1)*A(k-jmax,5) + psi*(A(k,2)*A(k-1,5)+ A(k,1)*A(k-jmax,4));
            A(k,4) = (A(k,4) - psi*A(k,1)*A(k-jmax,4))/A(k,3);
        end
        
        k=N;
        A(k,1) = A(k,1)/(1+psi*A(k-jmax,4));
        A(k,2) = A(k,2)/(1+psi*A(k-1,5));
        A(k,3) = A(k,3) - A(k,2)*A(k-1,4) - A(k,1)*A(k-jmax,5) + psi*(A(k,2)*A(k-1,5)+ A(k,1)*A(k-jmax,4));           
        
        for i=1:N
            for j=1:5
                LU(i,j) = A(i,j);
            end
        end
        
    case 9
        
        LU = zeros(N,stencil);
        
        for k = 1:N
            
            %a(k)
            if k>=jmax+2                
                LU(k,1) = A(k,1);
            end
            
            %b(k)
            if k>=jmax+1
                n = A(k,2);
                d = 1;
                if k-jmax-1>0
                    n = n - LU(k,1)*LU(k-jmax-1,6);
                end
                if k-jmax+1>0
                    n = n - psi*LU(k-jmax+1,6)*A(k,3);
                end
                if k-jmax>0
                    d = d - psi*LU(k-jmax,6)*LU(k-jmax+1,6);
                end
                LU(k,2) = n/d;
            end
            
            %c(k)
            if k>=jmax
               LU(k,3) = A(k,3);
               if k-jmax>0
                    LU(k,3) = LU(k,3) - LU(k,2)*LU(k-jmax,6);
               end
            end
            
            %d(k)
            if k>=2
                n = A(k,4);
                if k-jmax-1>0
                    n = n - 2*psi*LU(k,1)*LU(k-jmax-1,7);
                    n = n - LU(k,1)*LU(k-jmax-1,8);
                end
                if k-jmax>0
                    n = n - LU(k,2)*LU(k-jmax,7);
                end
                d = 1;
                if k-1>0
                    d = d + 2*psi*LU(k-1,7);
                end
                LU(k,4) = n/d;
            end
            
            %e(k)
            LU(k,5) = A(k,5);
            if k-jmax-1>0
                LU(k,5) = LU(k,5) - LU(k,1)*LU(k-jmax-1,9) + psi*LU(k,1)*LU(k-jmax-1,7);
            end
            if k-jmax>0
                LU(k,5) = LU(k,5) - LU(k,2)*LU(k-jmax,8);
            end
            if k-jmax+1>0
                LU(k,5) = LU(k,5) - LU(k,3)*LU(k-jmax+1,7) + 2*psi*LU(k,3)*LU(k-jmax+1,6) + psi*LU(k,3)*LU(k-jmax+1,9);
            end
            if k-1>0
                LU(k,5) = LU(k,5) - LU(k,4)*LU(k-1,6) + 2*psi*LU(k,4)*LU(k-1,7);
            end
            
            %f(k)
            if k<=N-1
                LU(k,6) = A(k,6);
                if k-jmax+1>0
                    LU(k,6) = LU(k,6) - 2*psi*LU(k,3)*(LU(k-jmax+1,6) + LU(k-jmax+1,9)) - LU(k,3)*LU(k-jmax+1,8);
                end
                if k-jmax>0
                    LU(k,6) = LU(k,6) - LU(k,2)*LU(k-jmax,9);
                end
                LU(k,6) = LU(k,6)/LU(k,5);
            end
            
            %g(k)
            if k<=N-jmax+1
                LU(k,7) = A(k,7);
                if k-1>0
                    LU(k,7) = LU(k,7) - LU(k,4)*LU(k-1,8);
                end
                LU(k,7) = LU(k,7)/LU(k,5);
            end
            
            %h(k)
            if k<=N-jmax
                LU(k,8) = A(k,8);
                if k-1>0
                    LU(k,8) = LU(k,8) - psi*LU(k,4)*LU(k-1,7) - LU(k,4)*LU(k-1,9);
                end
                LU(k,8) = LU(k,8)/LU(k,5);
            end
            
            if k<=N-jmax-1
                LU(k,9) = A(k,9)/LU(k,5);
            end
            
        end
         
end

end