function score = ADM(refImg, tesImg)
%========================================================================
%Copyright(c) 2011 Li Songnan
%All Rights Reserved.
%
%This is an implementation of the image quality metric described in
%the following paper:
%
%[1] S. Li, F. Zhang, L. Ma, K.N. Ngan, "Image Quality Assessment by Separately
%Evaluating Detail Losses and Additive Impairments", IEEE Trans. on Multimedia, accepted, 2011
%
%Kindly report any corrections to snli@ee.cuhk.edu.hk
%
%Input : (1) refImg: the reference image
%        (2) tesImg: the test image
%Output: (1) score:  the quality prediction
%========================================================================

% Translate color images to gray images
if ndims(refImg) == 3
    refImg = rgb2gray(refImg);
    tesImg = rgb2gray(tesImg);
end
refImg = double(refImg);
tesImg = double(tesImg);

% Wavelet transform
dwtmode('per','nodisp');
[O,Z] = wavedec2(refImg,4,'db2');
[T,Z] = wavedec2(tesImg,4,'db2');

% Decouple the restored image and
% the additive impairment image
[R,A] = decouple(O,T,Z);

% Contrast sensitivity function
viewDis = 4;
R2 = CSF(R,Z,viewDis);
A2 = CSF(A,Z,viewDis);
O2 = CSF(O,Z,viewDis);

% Contrast masking
MTR = get_CM_thresh(R2,Z);
MTA = get_CM_thresh(A2,Z);
R3 = CM(R2,MTA);
A3 = CM(A2,MTR);
O3 = abs(O2);

% Two quality measures: q1 and q2
q1=DLM(R3,O3,Z);
q2=AIM(A3,Z);

% Combination of the q1 and q2
a1=-0.815;
a2=1375;
score=q1+a1*(0.5-1/(1+exp(a2*q2)));
end

%=========================
% Decoupling additive impairments and detail losses
%=========================
function [R,A]=decouple(O,T,Z)
R=zeros(size(O));

% appromixation subband
lamda=4;theta=1;
[s,e] = get_subband_se(lamda,theta,Z);
R(s:e) = T(s:e);

% high frequency subbands
for lambda=1:4
    [sH,eH] = get_subband_se(lambda,2,Z);
    [sV,eV] = get_subband_se(lambda,3,Z);
    [sD,eD] = get_subband_se(lambda,4,Z);
    
    OH = O(sH:eH);
    TH = T(sH:eH);
    OV = O(sV:eV);
    TV = T(sV:eV);
    OD = O(sD:eD);
    TD = T(sD:eD);
    
    kH=min(max(TH./(OH+1e-30),0),1);
    kV=min(max(TV./(OV+1e-30),0),1);
    kD=min(max(TD./(OD+1e-30),0),1);
    tmpH = kH.*OH;
    tmpV = kV.*OV;
    tmpD = kD.*OD;
    
    % special case: contrast enhancement
    OA = cal_angle(OH,OV);
    TA = cal_angle(TH,TV);
    diff = abs(OA-TA)*180/pi;
    pos = find(diff<1);
    tmpH(pos) = TH(pos);
    tmpV(pos) = TV(pos);
    tmpD(pos) = TD(pos);
    
    R(sH:eH) = tmpH;
    R(sV:eV) = tmpV;
    R(sD:eD) = tmpD;
end

% calcute DWT coeff. of 
% the additive impairment image
A=T-R;
end

%=========================
% Get the DWT coeff. indices of the DWT subband {lambda, theta} 
%=========================
function [s,e] = get_subband_se(lambda,theta,Z)
if lambda==4 && theta==1
    s=1; e=Z(1,1)*Z(1,2);
    return;
end
lambda2 = 5-lambda;
s=Z(1,1)*Z(1,2)+1;
for i=2:lambda2
    s=s+3*Z(i,1)*Z(i,2);
end
theta2 = theta-2;
s=s+theta2*Z(lambda2+1,1)*Z(lambda2+1,2);
e=s+Z(lambda2+1,1)*Z(lambda2+1,2)-1;
end

%=========================
% Implementation of Eq. (13) of [1]
%=========================
function OA=cal_angle(OH,OV)
OA = atan(OV./(OH+1e-30))+pi*(OH<0);
end

%=========================
% Cosntrast sensitivity function
%=========================
function X2 = CSF(X,Z,viewDis)
X2=zeros(size(X));
fs=Z(6,1);
[s,e] = get_subband_se(4,1,Z);
X2(s:e)=H(0)*X(s:e);
for lambda=1:4
    f1=pi*fs*viewDis/(180*2^lambda);
    for theta=2:4
        p = [1,1,-1];
        f2 = f1/(0.15*p(theta-1)+0.85);
        [s,e] = get_subband_se(lambda,theta,Z);
        X2(s:e) = H(f2).*X(s:e);
    end
end
end

%=========================
% Implementation of Eq. (16) of [1]
%=========================
function c=H(w)
a=0.31;
b=0.69;
c=0.29;
c=(a+b*w)*exp(-c*w);
end

%=========================
% Calculate the treshold values for contrast masking
%=========================
function MT = get_CM_thresh(M,Z)
MT = zeros(size(M));

% high frequency subbands only
for lambda=1:4
    [sH,eH] = get_subband_se(lambda,2,Z);
    [sD,eD] = get_subband_se(lambda,4,Z);
    
    H = get_subband_2D(M,lambda,2,Z);
    V = get_subband_2D(M,lambda,3,Z);
    D = get_subband_2D(M,lambda,4,Z);
    
    w=[1/30,1/30,1/30;1/30,1/15,1/30;1/30,1/30,1/30];
    MH = conv2(abs(H),w,'same');
    MV = conv2(abs(V),w,'same');
    MD = conv2(abs(D),w,'same');
    
    MTlambda=MH(:)+MV(:)+MD(:);
    MT(sH:eD) = [MTlambda;MTlambda;MTlambda];
end
end

%=========================
% Get the DWT coeff. indices of certain DWT subband
% and then reshape it into a 2D map
%=========================
function Xsubband2D = get_subband_2D(X,lambda,theta,Z)
h = Z(6-lambda,1);
w = Z(6-lambda,2);
[s,e] = get_subband_se(lambda,theta,Z);
Xsubband2D=reshape(X(s:e),[h,w]);
end

%=========================
% Implementation of the contrast masking
%=========================
function X3 = CM(X2,MT)
X2=abs(X2);
X3=max((X2-MT),0);
end

%=========================
% The detail loss measure
%=========================
function q1=DLM(R3,O3,Z)
nume=0;
deno=0;
for lambda=1:4
    for theta=2:4
        Rsub = get_subband_2D(R3,lambda,theta,Z);
        Osub = get_subband_2D(O3,lambda,theta,Z);
        RsubCen = center(Rsub);
        OsubCen = center(Osub);
        a=RsubCen.^3;
        a=sum(a(:))^(1/3);
        b=OsubCen.^3;
        b=sum(b(:))^(1/3);
        nume=nume+a;
        deno=deno+b;
    end
end

q1=nume/deno;
end

%=========================
% The additive impairment measure
%=========================
function q2=AIM(A3,Z)
nume=0;
deno=numel(A3);
for lambda=1:4
    for theta=2:4
        Asub = get_subband_2D(A3,lambda,theta,Z);
        AsubCen = center(Asub);
        a=AsubCen.^3;
        a=sum(a(:))^(1/3);
        nume=nume+a;
    end
end

q2=nume/deno;
end

%=========================
% Take the center part of an image
%=========================
function B=center(A)
marg=0.1;
[h,w]=size(A);
margH = round(max(1,h*marg));
margW = round(max(1,w*marg));
B = A(margH:h-margH+1,margW:w-margW+1);
end
