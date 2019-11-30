%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% unit vector between two points in caertesian frame
%
%%% inputs
% ----------
% X : array, size(3,1), point in cartesian frame
% Y : array, size(3,1), point in cartesian frame
%
%%% outputs
%-----------
% e : array, size(3,1), unit vector from X to Y
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function e = unit_vector(X,Y)
    % unit vector from X to Y
    e = (Y - X)./norm(Y - X);
end

