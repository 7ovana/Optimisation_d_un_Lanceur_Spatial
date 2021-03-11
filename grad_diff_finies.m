function [gradf, gradc] = grad_diff_finies(x, h, problem)
     [f, c] = problem(x); 
     n = length(x);
     m = length(c);
     % on initialise grad et jac_c: 
     gradf = zeros(n,1);
     gradc = zeros(m,n);
     %Matrice des déplacements h:
     H = diag(h);
  
     for i = 1:n
       %différences finies, version décentrée:
       [fh, ch] = problem(x + H(:,i));
       gradf(i) = (fh - f) / h(i);
       gradc(:,i) = (ch - c) ./ h(i);  
     end
     
     gradc = gradc';
end
