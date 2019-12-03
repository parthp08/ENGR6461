%%% calculation for 1-sigma error
%%% have to compte the DOP matrix in ENU frame
%%% Note: have to run main.m file before running this file
%
%%% References
% -------------
%[1] 'Global Positioning System: Theory and Applications'; edited by 
%     Parkinson, Spilker, Axelrad, Enge; AIAA; 1996
%[2] 'Aerospace Navigation Systems'; edited by A.V.Nebylv, J.Watson
%[3] https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
%[4] https://stattrek.com/online-calculator/chi-square.aspx

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

UERE = 5.3; % page 481 Ref[1] % 1-sigma user equivalent error

% for east and north plane error
ENDOP = DOP_enu(1:2,1:2);

% covarience matrix
K = (UERE^2)*ENDOP; % it is correlated  %%%%%%%% put UERE^2
sigma_xe = sqrt(K(1,1));
sigma_xn = sqrt(K(2,2));

S = 2.28; % scale of the ellipse for 1-sigma (68.2% confidence)
[V,E] = eig(ENDOP);
a = sigma_xe*sqrt(S); % semi-major axis % in north direction
b = sigma_xn*sqrt(S); % semi-minor axis % in east direction
theta = atan(V(2,2)/V(2,1)); % orientation of the ellipse

xe0=0; % x0,y0 ellipse centre coordinates
xn0=0;
t=-pi:0.01:pi;
XE=xe0+a*cos(t);
XN=xn0+b*sin(t);
% rotate ellipse to theta angle
XE_rot = XE*cos(theta) + -XN*sin(theta);
XN_rot = XE*(sin(theta)) + XN*cos(theta);
plot(XE_rot,XN_rot);
title("1-sigma ellipse");
xlabel("X east");
ylabel("X north");
axis equal;
