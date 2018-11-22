% Li Songnan, XX-XX-2014, The Chinese University of Hong Kong
%
% This is an implementation of the following paper:
% [1] Songnan Li, King Ngi Ngan, Lu Sheng, "Screen-camera calibration using a
% thread", submitted to ICIP 2014
% 
% For measuring color differences between a pixel and its neighbors, we use this code:
% http://www.mathworks.com/matlabcentral/answers/73741#comment_145951
%
% For the RANSAC algorithm, we use this toolbox:
% https://github.com/RANSAC/RANSAC-Toolbox
%
% Please Kindly report any corrections to snli@ee.cuhk.edu.hk

function screen_camera_cali
clc;
clear;

% include the RANSAC Toolbox
addpath('./Ransac/');

% load the default values (inputs of the last run) 
load tmp.mat

% camera intrisic parameters
% no input then use the default value
FXInput = input('focal length in X direction: ');
if ~isempty(FXInput) FX = FXInput; end;   
disp(sprintf('%.2f',FX));

FYInput = input('focal length in Y direction: ');
if ~isempty(FYInput) FY = FYInput; end;
disp(sprintf('%.2f',FY));

TXInput = input('image center of X direction: ');
if ~isempty(TXInput) TX = TXInput; end;
disp(sprintf('%.2f',TX));

TYInput = input('image center of Y direction: ');
if ~isempty(TYInput) TY = TYInput; end;
disp(sprintf('%.2f',TY));

% we assume that the skew coefficient is zero
IntrinsicMatrix = [FX,0,TX;0,FY,TY;0,0,1];

% width and height of the screen
widthScreenInput = input('Width of the screen (cm): ');
if ~isempty(widthScreenInput)  widthScreen  = widthScreenInput; end;
disp(sprintf('%.2f',widthScreen));

heightScreenInput = input('Height of the screen (cm): ');
if ~isempty(heightScreenInput) heightScreen = heightScreenInput; end;
disp(sprintf('%.2f',heightScreen));

% basename of the calibration images
basenameInput = input('Basename calibration images (without number nor suffix): ','s');
if ~isempty(basenameInput)  basename  = basenameInput; end;
disp(basename);

% save the user's inputs
save tmp.mat FX FY TX TY widthScreen heightScreen basename

% upper-left screen corner
alphaBeta1 = detect_lines_4_each_corner('./1_upperleft/',basename);
show_lines(alphaBeta1,-20000,20000,'--y');                             
% reject outliers using Ransac
[alphaBeta1 m1] = reject_outliers(alphaBeta1);
show_lines(alphaBeta1,-20000,20000,'-r');                       
% show the intersection point
figure(2); hold on; xlabel('x'); ylabel('y');
plot(m1(1),m1(2),'*c');pause(1);

% upper-right screen corner
alphaBeta2 = detect_lines_4_each_corner('./2_upperright/',basename);
show_lines(alphaBeta2,-20000,20000,'--y');                           
% reject outliers using Ransac
[alphaBeta2 m2] = reject_outliers(alphaBeta2);
show_lines(alphaBeta2,-20000,20000,'-g');                     
% show the intersection point
plot(m2(1),m2(2),'*c');pause(1);

% bottom-left screen corner
alphaBeta3 = detect_lines_4_each_corner('./3_bottomleft/',basename);
show_lines(alphaBeta3,-20000,20000,'--y');                         
% reject outliers using Ransac
[alphaBeta3 m3] = reject_outliers(alphaBeta3);
show_lines(alphaBeta3,-20000,20000,'-b');                   
% show the intersection point
plot(m3(1),m3(2),'*c');pause(1);

% bottom-right screen corner
alphaBeta4 = detect_lines_4_each_corner('./4_bottomright/',basename);
show_lines(alphaBeta4,-20000,20000,'--y');                        
% reject outliers using Ransac
[alphaBeta4 m4] = reject_outliers(alphaBeta4);
show_lines(alphaBeta4,-20000,20000,'-m');                 
% show the intersection point
hold on; plot(m4(1),m4(2),'*c');pause(1);

% initialize the screen position
w = widthScreen;
h = heightScreen;
m = [m1,m2,m3,m4];
M1 = [0;0];
M2 = [w;0];
M3 = [0;h];
M4 = [w;h];
M = [M1,M2,M3,M4];
[P1 thetaX thetaY thetaZ] = ini_par(m,M,IntrinsicMatrix);

% minimize projection errors
opt = optimset('algorithm', 'levenberg-marquardt',...       % Use L-M algorithm for optimization
    'tolfun', 1e-8, ...
    'tolx', 1e-8, ...
    'maxfunevals', 100000, ...
    'MaxIter', 10000, ...
    'display', 'off');

XandRinit = [P1;thetaX;thetaY;thetaZ];
XandR = lsqnonlin(@(x) projection_err(x,h,w,IntrinsicMatrix,alphaBeta1,alphaBeta2,alphaBeta3,alphaBeta4),...
    XandRinit, [], [], opt);
P1 = XandR(1:3);
R = get_rotation(XandR(4),XandR(5),XandR(6));
P2 = P1+w*R*[1;0;0];
P3 = P1+h*R*[0;1;0];
P4 = P1+w*R*[1;0;0]+h*R*[0;1;0];

% draw projections on the image plane
proj1 = get_projection(IntrinsicMatrix,P1);
proj2 = get_projection(IntrinsicMatrix,P2);
proj3 = get_projection(IntrinsicMatrix,P3);
proj4 = get_projection(IntrinsicMatrix,P4);
figure(2);
plot(proj1(1),proj1(2),'oc');
plot(proj2(1),proj2(2),'oc');
plot(proj3(1),proj3(2),'oc');
plot(proj4(1),proj4(2),'oc');
pause(1);

% draw screen position
figure(3);hold on;grid on;
xlabel('x');ylabel('y');zlabel('z');
plot3([P1(1);P2(1);P4(1);P3(1);P1(1)],...
      [P1(2);P2(2);P4(2);P3(2);P1(2)],...
      [P1(3);P2(3);P4(3);P3(3);P1(3)],'k');
C = (P1+P2+P3+P4)/4;
axis([-40+C(1), 40+C(1), -40+C(2), 40+C(2), -10, 70]);
view([120,30]);
text(P1(1),P1(2),P1(3),['(' sprintf('%.2f',P1(1)) ',' sprintf('%.2f',P1(2)) ',' sprintf('%.2f',P1(3)) ')']);
text(P2(1),P2(2),P2(3),['(' sprintf('%.2f',P2(1)) ',' sprintf('%.2f',P2(2)) ',' sprintf('%.2f',P2(3)) ')']);
text(P3(1),P3(2),P3(3),['(' sprintf('%.2f',P3(1)) ',' sprintf('%.2f',P3(2)) ',' sprintf('%.2f',P3(3)) ')']);
text(P4(1),P4(2),P4(3),['(' sprintf('%.2f',P4(1)) ',' sprintf('%.2f',P4(2)) ',' sprintf('%.2f',P4(3)) ')']);
P1 = P1'
P2 = P2'
P3 = P3'
P4 = P4'

% draw XYZ axis
plot3([0,4],[0,0],[0,0],'-r','LineWidth',3);
plot3([0,0],[0,4],[0,0],'-g','LineWidth',3);
plot3([0,0],[0,0],[0,4],'-b','LineWidth',3);

% draw camera position
camra1=[-4;4;4];
camra2=[4;4;4];
camra3=[-4;-4;4];
camra4=[4;-4;4];
plot3([camra1(1);camra2(1);camra4(1);camra3(1);camra1(1)],...
      [camra1(2);camra2(2);camra4(2);camra3(2);camra1(2)],...
      [camra1(3);camra2(3);camra4(3);camra3(3);camra1(3)],'m');
plot3([0,camra1(1)],[0,camra1(2)],[0,camra1(3)],'-m','LineWidth',2);
plot3([0,camra2(1)],[0,camra2(2)],[0,camra2(3)],'-m','LineWidth',2);
plot3([0,camra3(1)],[0,camra3(2)],[0,camra3(3)],'-m','LineWidth',2);
plot3([0,camra4(1)],[0,camra4(2)],[0,camra4(3)],'-m','LineWidth',2);

end

% P2 and P1 has this relationship: P2 = P1 + w*R*[1,0,0]'
% P2 should be as close to "Ray2" as possible
% P2 has a closed form solution as given below
function P2 = find_P2(P1,Ray2)
D = [-1;0;0];
coeff = inv([D'*D,-Ray2*D; D'*Ray2',-Ray2*Ray2']) * [-P1'*D;-P1'*Ray2'];
P2 = P1+coeff(1)*D;
end

% Transform the intersection point into a normalized 3D ray
function ray3D = trans_point_to_3DRay(m,IntrinsicMatrix)
row1 = IntrinsicMatrix(1,:);
row2 = IntrinsicMatrix(2,:);
row3 = IntrinsicMatrix(3,:);
v1 = m(1)*row3-row1;
v2 = m(2)*row3-row2;
ray3Dtmp = cross(v1,v2);
ray3D = ray3Dtmp/norm(ray3Dtmp);
end

% Reject outliers using RANSAC
function  [alphaBetaArr2, m] = reject_outliers(alphaBetaArr)
% Set RANSAC options
sigma = 15;
options.epsilon = 1e-6;
options.P_inlier = 0.99;
options.sigma = sigma;
options.est_fun = @estimate_m;
options.man_fun = @error_m;
options.mode = 'MSAC';
options.Ps = [];
options.notify_iters = [];
options.min_iters = 1000;
options.fix_seed = false;
options.reestimate = true;
options.stabilize = false;

% Perform RANSAC to get inliers
X = alphaBetaArr';
results = RANSAC(X, options);

% Intersection point
m = results.Theta;

% Reject outliers
inliersIdx = results.CS;
alphaBetaArr2 = [];
for i=1:size(alphaBetaArr,1)
    if(inliersIdx(i)==1)
        alphaBetaArr2 = [alphaBetaArr2;alphaBetaArr(i,:)];
    end
end

end

% Calculate the intersection point of multiple lines
function [Theta, k] = estimate_m(X, s)
k = 2; % at least 2 lines are required

if (nargin == 0) || isempty(X)
    Theta = [];
    return;
end;

if (nargin == 2) && ~isempty(s)
    X = X(:, s);
end;

% check if we have enough points
N = size(X, 2);
if (N < k)
    error('estimate_m:inputError');
end;

sumAiAit=zeros(2,2);
sumciAi=zeros(2,1);
for pp = 1:size(X,2)
    alpha = X(1,pp);
    beta = X(2,pp);
    
    Ai = [alpha;-1];
    ci = beta;
    sumAiAit=sumAiAit+Ai*Ai';
    sumciAi=sumciAi+ci*Ai;
end
Theta = inv(sumAiAit + eye(2)*1e-10)*(-sumciAi); % add small values to avoid matrix singularity

return;
end

function [E T_noise_squared d] = error_m(Theta, X, sigma, P_inlier)
% compute the squared error
E = [];
if ~isempty(Theta) && ~isempty(X)
    E = ( Theta(1)*X(1,:) - Theta(2) + X(2,:) ).^2 ./ (X(1,:).^2 + 1);        
end;

% compute the error threshold
if (nargout > 1)
    
    if (P_inlier == 0)
        T_noise_squared = sigma;
    else
        % Assumes the errors are normally distributed. Hence the sum of
        % their squares is Chi distributed (with 2 DOF since we are 
        % computing the distance of a 2D point to a line)
        d = 2;
        
        % compute the inverse probability
        T_noise_squared = sigma^2 * chi2inv_LUT(P_inlier, d);

    end;
    
end;

end

% Show the defected lines on the image plane
function show_lines(alphaBetaArr,minX,maxX,style)
for i = 1:size(alphaBetaArr,1)
        x=minX:maxX;
        alpha = alphaBetaArr(i,1);
        beta = alphaBetaArr(i,2);
        y=alpha*x+beta;
        figure(2);hold on;plot(x,y,style);
end
end

function alphaBetaArr = detect_lines_4_each_corner(subfolder,basename)
l = dir([subfolder basename '*']);
Nl = size(l,1);
alphaBetaArr = [];
pp = 1;
manuallySetFlag = false;
while pp<=Nl
        fileName =  [subfolder l(pp).name];
        [alpha, beta, noLineDetectedFlag] = get_line(fileName,manuallySetFlag);
        if noLineDetectedFlag == true; continue; end;
        manuallySetFlag = false;
        
        key = input('[]=Next, 1=Redo,  2=Manually,  3=Skip:');
        if isempty(key)
            % use the line detection result (default)
            alphaBetaArr = [alphaBetaArr;alpha, beta];
        elseif key == 1
            % redo the line detection
            pp = pp - 1;
        elseif key == 2;
            % mannually set the line position
            manuallySetFlag = true;
            pp = pp - 1;
        elseif key == 3;
            % skip this frame
            disp('Frame is skipped');
        elseif noLineDetectedFlag == true
            pp = pp - 1;
        end
        
pp = pp+1;
end

end

% Semi-automatic line detection
function [alpha, beta, noLineDetectedFlag] = get_line(fileName,manuallySetFlag)
% read image and transform to LAB color space
rgbImage  = imread(fileName);
cform = makecform('srgb2lab');
lab_Image = applycform(im2double(rgbImage),cform);

% extract out the color bands from the original image
% into 3 separate 2D arrays, one for each color component.
LChannel = lab_Image(:, :, 1);
aChannel = lab_Image(:, :, 2);
bChannel = lab_Image(:, :, 3);
windowSize = 8;
[LMean, aMean, bMean] = GetMeanLABValues(LChannel, aChannel, bChannel, windowSize);

% get the elta colors for each LAB color channel.
deltaL = LMean - LChannel;
deltaa = aMean - aChannel;
deltab = bMean - bChannel;

% create the Delta E image.
% this is an image that represents the color difference.
% Delta E is the square root of the sum of the squares of the delta images.
deltaE = sqrt(deltaL .^ 2 + deltaa .^ 2 + deltab .^ 2);

% thresholding
BW =logical(zeros(size(deltaE)));
BW(find(deltaE>8)) = 1;

% rough line position given by the user
figure(1);hold off;imshow(rgbImage);
[x,y] = ginput(1);
x1=x(1);y1=y(1);
figure(1);hold on;
plot(x1,y1,'+','color',[ 1.000 0.314 0.510 ],'linewidth',2);
[x,y] = ginput(1);
x2=x(1);y2=y(1);
plot(x2,y2,'+','color',[ 1.000 0.314 0.510 ],'linewidth',2);

% get alpha and beta (y = alpha*x + beta)
alpha = (y2-y1)/(x2-x1);
beta = -(y2-y1)*x1/(x2-x1)+y1;

if manuallySetFlag == true
    % show manually determined line
    x = 1:size(rgbImage,2);
    y =  alpha*x + beta;
    figure(1); plot(x,y,'-g');
    
    noLineDetectedFlag = false;
    return;
end

% get theta and rho
theta = acot(-alpha)*180/pi;
rho = beta*sin(theta*pi/180);
theta = round(theta);    % round to the nearest integer
rho = round(rho);        % round to the nearest integer

% Hough transform
[H,T,R] = hough(BW);

% refine the line position
maxEdgeNo = 0;
searchRangeT = 6;
searchRangeR = 6;
pT = find(T==theta);
pR = find(R==rho);
for i=-searchRangeR:searchRangeR
    for j=-searchRangeT:searchRangeT
        rhotmp = pR+i;
        thetatmp = pT+j;
        if rhotmp > 0 && rhotmp <= size(H,1) && thetatmp > 0 && thetatmp <= size(H,2) && H(rhotmp,thetatmp)~=0
            if (H(rhotmp,thetatmp)>maxEdgeNo)
                maxEdgeNo = H(rhotmp,thetatmp);
                P = [rhotmp,thetatmp];
            end
        end
    end
end

% Hough line detection
lines = houghlines(BW,T,R,P,'FillGap',5,'MinLength',7);
if isfield(lines,'rho') == false
    noLineDetectedFlag = true;
    disp('No line detected');
    return;
else
    noLineDetectedFlag = false;
end
rhoDetected = lines(1).rho;
thetaDetected = lines(1).theta;
alpha = -cot(thetaDetected*pi/180);
beta = rhoDetected/sin(thetaDetected*pi/180);

% show detected line
x = 1:size(rgbImage,2);
y =  alpha*x + beta;
figure(1); plot(x,y,'-g');

end

% Get the mean LAB values
function [LMean, aMean, bMean] = GetMeanLABValues(LChannel, aChannel, bChannel, windowSize)
try
	kernel = ones(windowSize, windowSize) / windowSize^2;
	LMean = conv2(LChannel, kernel, 'same'); % Average of only the pixels within the windowed area.
	aMean = conv2(aChannel, kernel, 'same'); % Average of only the pixels within the windowed area.
	bMean = conv2(bChannel, kernel, 'same'); % Average of only the pixels within the windowed area.
catch ME
	errorMessage = sprintf('Error running GetMeanLABValues:\n\n\nThe error message is:\n%s', ...
		ME.message);
	WarnUser(errorMessage);
end

end

% Get the rotation matrix
function R = get_rotation(thetaX,thetaY,thetaZ)
rotX=[1,              0,                0;...
      0,              cos(thetaX),      -sin(thetaX);...
      0,              sin(thetaX),      cos(thetaX)];
rotY=[cos(thetaY),    0,                sin(thetaY);...
      0,              1,                0;...
      -sin(thetaY),   0,                cos(thetaY)];
rotZ=[cos(thetaZ),    -sin(thetaZ),     0;...
      sin(thetaZ),    cos(thetaZ),      0;...
      0,              0,                1];
R=rotZ*rotY*rotX;
end

% Decompose the rotation matrix
function [x,y,z] = decompose_rotation(R)
	x = atan2(R(3,2), R(3,3));
	y = atan2(-R(3,1), sqrt(R(3,2)*R(3,2) + R(3,3)*R(3,3)));
	z = atan2(R(2,1), R(1,1));
end

% Calculate the projection errors
function cost = projection_err(XandR,H,W,intriMax,alphaBeta1,alphaBeta2,alphaBeta3,alphaBeta4)
cost = [];
P1 = XandR(1:3);
R = get_rotation(XandR(4),XandR(5),XandR(6));
P2 = P1+W*R*[1;0;0];
P3 = P1+H*R*[0;1;0];
P4 = P1+W*R*[1;0;0]+H*R*[0;1;0];
pro1 = get_projection(intriMax,P1);
pro2 = get_projection(intriMax,P2);
pro3 = get_projection(intriMax,P3);
pro4 = get_projection(intriMax,P4);

for i=1:size(alphaBeta1,1)
    alpha=alphaBeta1(i,1);
    beta=alphaBeta1(i,2);

    a=alpha;
    b=-1;
    c=beta;
    
    dis = distance2D_Point2Plane(pro1,a,b,c);
    cost = [cost; dis/sqrt(size(alphaBeta1,1))];
end

for i=1:size(alphaBeta2,1)
    alpha=alphaBeta2(i,1);
    beta=alphaBeta2(i,2);
   
    a=alpha;
    b=-1;
    c=beta;
    
    dis = distance2D_Point2Plane(pro2,a,b,c);
    cost = [cost; dis/sqrt(size(alphaBeta2,1))];
end

for i=1:size(alphaBeta3,1)
    alpha=alphaBeta3(i,1);
    beta=alphaBeta3(i,2);
    
    a=alpha;
    b=-1;
    c=beta;
    
    dis = distance2D_Point2Plane(pro3,a,b,c);
    cost = [cost; dis/sqrt(size(alphaBeta3,1))];
end

for i=1:size(alphaBeta4,1)
    alpha=alphaBeta4(i,1);
    beta=alphaBeta4(i,2);
   
    a=alpha;
    b=-1;
    c=beta;
    
    dis = distance2D_Point2Plane(pro4,a,b,c);
    cost = [cost; dis/sqrt(size(alphaBeta4,1))];
end

end

% Calculate the projection errors (enforcing no roll rotation)
function cost = projection_err_zero_roll(XandR,H,W,intriMax,alphaBeta1,alphaBeta2,alphaBeta3,alphaBeta4)
cost = [];
P1 = XandR(1:3);
R = get_rotation(XandR(4),XandR(5),0);
P2 = P1+W*R*[1;0;0];
P3 = P1+H*R*[0;1;0];
P4 = P1+W*R*[1;0;0]+H*R*[0;1;0];
pro1 = get_projection(intriMax,P1);
pro2 = get_projection(intriMax,P2);
pro3 = get_projection(intriMax,P3);
pro4 = get_projection(intriMax,P4);

for i=1:size(alphaBeta1,1)
    alpha=alphaBeta1(i,1);
    beta=alphaBeta1(i,2);

    a=alpha;
    b=-1;
    c=beta;
    
    dis = distance2D_Point2Plane(pro1,a,b,c);
    cost = [cost; dis/sqrt(size(alphaBeta1,1))];
end

for i=1:size(alphaBeta2,1)
    alpha=alphaBeta2(i,1);
    beta=alphaBeta2(i,2);
   
    a=alpha;
    b=-1;
    c=beta;
    
    dis = distance2D_Point2Plane(pro2,a,b,c);
    cost = [cost; dis/sqrt(size(alphaBeta2,1))];
end

for i=1:size(alphaBeta3,1)
    alpha=alphaBeta3(i,1);
    beta=alphaBeta3(i,2);
    
    a=alpha;
    b=-1;
    c=beta;
    
    dis = distance2D_Point2Plane(pro3,a,b,c);
    cost = [cost; dis/sqrt(size(alphaBeta3,1))];
end

for i=1:size(alphaBeta4,1)
    alpha=alphaBeta4(i,1);
    beta=alphaBeta4(i,2);
   
    a=alpha;
    b=-1;
    c=beta;
    
    dis = distance2D_Point2Plane(pro4,a,b,c);
    cost = [cost; dis/sqrt(size(alphaBeta4,1))];
end

end

% Calculate the perspective projection of a 3D point X
function m = get_projection(intriMax,X)
m=zeros(2,1);
Mtmp = intriMax*X;
m(1) = Mtmp(1)/Mtmp(3);
m(2) = Mtmp(2)/Mtmp(3);
end

% Calculate the distance between a 2D point and a 2D line
function D = distance2D_Point2Plane(x,a,b,c)
D = abs(a*x(1)+b*x(2)+c)/sqrt(a^2+b^2);
end

function [P1 thetaX thetaY thetaZ] = ini_par(m,M,A)
% Compute homography matrix 
H = compute_homography(m, M);
h1 = H(:,1);
h2 = H(:,2);
h3 = H(:,3);

% Calculate the camera's roation and translation matrix
lambda = 1/norm(inv(A)*h1);
r1 = lambda*inv(A)*h1;
r2 = lambda*inv(A)*h2;
r3 = cross(r1,r2);
t = lambda*inv(A)*h3;

% Refine the rotation matrix 
R = [r1,r2,r3];
[U,S,V] = svd(R);
R = U*V';
[x,y,z] = decompose_rotation(R);

P1 = t;
thetaX=x;
thetaY=y;
thetaZ=z;
end

function [H,Hnorm,inv_Hnorm] = compute_homography(m, M)
% compute the planar homography between the point coordinates on the plane
% (M) and the point coordinates (m)
Np = size(m, 2);
% complete matrices
if size(m, 1) < 3
    m = [m; ones(1, Np)];
end
if size(M, 1) < 3
    M = [M; ones(1, Np)];
end

m = m./(ones(3,1)*m(3,:));
M = M./(ones(3,1)*M(3,:));

% pre-normalized , important
ax = m(1,:); ay = m(2,:);
mx = mean(ax); my = mean(ay);

ax = ax - mx; ay = ay - my;
scx = mean(abs(ax)); scy = mean(abs(ay));

Hnorm = [1/scx, 0, -mx/scx; 0, 1/scy, -my/scy; 0, 0, 1];
inv_Hnorm = [scx, 0, mx; 0, scy, my; 0, 0, 1];
mn = Hnorm*m;

% compute the homography
L = zeros(2*Np, 9);
L(1:2:2*Np, 1:3) = M';
L(1:2:2*Np, 7:9) = -((ones(3,1)*mn(1,:)).*M)';  
L(2:2:2*Np, 4:6) = M';
L(2:2:2*Np, 7:9) = -((ones(3,1)*mn(2,:)).*M)';

if Np > 4
    L = L'*L;
end

[U,S,V] = svd(L);
hh = V(:,9);
hh = hh/hh(9);
Hrem = reshape(hh, 3, 3)';

H = inv_Hnorm*Hrem;

end

