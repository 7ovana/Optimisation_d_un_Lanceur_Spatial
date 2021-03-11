function [res, nb_iter] = newton(x0, f, eps)
  
  x = x0;
  [y, f_prim] = f(x0);
  nb_iter = 0;
  
  while 1
    y_avant = y;
    x = x - y / f_prim;
    if (x < 1) 
        x = 1.1; 
    end
    [y, f_prim] = f(x);
    nb_iter = nb_iter + 1;  
    if (abs(y - y_avant) < eps) 
      break; 
    end
  end
  res = x;
end
