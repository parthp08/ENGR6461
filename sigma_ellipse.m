%%% calculation for 1-sigma error
%%% have to compte the DOP matrix in ENU frame
%%% 1-sigma ellipse means 68.2689492% of the time values will be 
%%% within the boundary of ellipse [5].
%
%%% References
% [1] https://www.mathworks.com/matlabcentral/answers/86615-how-to-plot-an-ellipse#
% [2] 'GNSS Error Sources',Malek Karaim, Mohamed Elsheikh and Aboelmagd Noureldin 
% [3] https://www.xarg.org/2018/04/how-to-plot-a-covariance-error-ellipse/
% [4] https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
% [5] https://en.wikipedia.org/wiki/Standard_deviation

% use Pu and P_sat_arr from main files and convert them to ENU Frame
Pu_enu = ECEF2ENU(Pu);
P_sat_arr_enu = ECEF2ENU(P_sat_arr);

%%%%% construct H matrix
% number of satellite
n = size(P_sat_arr_enu, 2); % 2 for size of columns
% vector dimension
p = size(P_sat_arr_enu(:,1), 1); % 1 for size of rows
% initialize H matrix
H = ones(n, p+1);
% fill H matrix
for i = 1:n
    H(i,1:p) = (-unit_vector(Pu_enu, P_sat_arr_enu(:,i)));
end

%%%%%% calculate DOP in ENU frame
H_T = H';
DOP_enu = inv(H_T*H);
GDOP_enu = sqrt(DOP_enu(1,1)+DOP_enu(2,2)+DOP_enu(3,3)+DOP_enu(4,4)); % geometric DOP
PDOP_enu = sqrt(DOP_enu(1,1)+DOP_enu(2,2)+DOP_enu(3,3)); % Position DOP
HDOP_enu = sqrt(DOP_enu(1,1)+DOP_enu(2,2)); % Horizontal DOP

%%%%%% 1-sigma ellipse  % Ref [1] and [3]
% 1-sigma ellipse in East and North Direction
x_e = DOP_enu(1,1); % semimajor axis of an ellipse
x_n = DOP_enu(2,2); % semiminor axis of an ellipse
xe0=0; % x0,y0 ellipse centre coordinates
xn0=0;
t=-pi:0.01:pi;
XE=xe0+x_e*cos(t);
XN=xn0+x_n*sin(t);
plot(XE,XN);
title("1-sigma ellipse");
xlabel("X east");
ylabel("X north");
axis equal;