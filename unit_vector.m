function e = unit_vector(X,Y)
    % unit vector from X to Y
   
    e = (Y - X)./norm(Y - X);
    
end

