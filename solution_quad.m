function [d,lambda]= solution_quad(gradf,gradc,c,Q)
 g = gradf;
 A = gradc';
 b = -c ;

 % lambda= -inv(A*inv(Q)* A')*(A * inv(Q)*g+b);
 %dk = - inv(Q)*(A' * lambda + g) ; 
 
 lambda = -(A*(Q\A'))\(A*(Q\g)+b);
 d = - Q\(A' * lambda + g) ;                    
end 