%=========================================================================
%This is an implementation of the video quality metric described in
%the following paper:
%
%   Songnan Li, Lin Ma, King Ngi Ngan, "Video quality assessment by
%   decoupling additive impairments and detail losses", The third International 
%   Workshop on Quality of Multimedia Experience, pp. 90-95, Sep. 2011.
%=========================================================================

function score = VADM(refAddr,disAddr,framenum,height,width,viewingDis,chrotype)
% refAddr:    address of the reference video
% disAddr:    address of the distorted video
% framenum:   total frame number
% height:     frame height 
% width:      frame width 
% viewingDis: viewing distance
% chrotype:   chroma type (handle "420" only)

framescoresDLM = zeros(1,framenum);
framescoresAIM = zeros(1,framenum);
for frIdx=1:framenum
    
    % Read current frame
    refFrame = get_one_frame(refAddr,frIdx,chrotype,height,width);
    disFrame = get_one_frame(disAddr,frIdx,chrotype,height,width);
    refFrameY_xn = refFrame{1}; % use Y component only
    disFrameY_xn = disFrame{1};
    if frIdx==1
        refFrameY_xn_minus_1 = refFrameY_xn;
        disFrameY_xn_minus_1 = disFrameY_xn;
        refFrameY_yn_minus_1 = refFrameY_xn;
        disFrameY_yn_minus_1 = disFrameY_xn;
    end
    
    % Temporal filtering
    refFrameY_yn=0.8*refFrameY_xn+0.12*refFrameY_xn_minus_1+0.08*refFrameY_yn_minus_1;
    disFrameY_yn=0.8*disFrameY_xn+0.12*disFrameY_xn_minus_1+0.08*disFrameY_yn_minus_1;
    refFrameY=refFrameY_yn;
    disFrameY=disFrameY_yn;
    
    % Decoupling
    dwtmode('per','nodisp');
    wavelet = 'db1';
    totalLevel = 4;
    refFrameY = crop_frame(refFrameY,32);
    disFrameY = crop_frame(disFrameY,32);
    [O,Z] = wavedec2(refFrameY,totalLevel,wavelet);
    [D,Z] = wavedec2(disFrameY,totalLevel,wavelet);
    [R,A]=decouple4video(O,D,Z,totalLevel);
    
    % CSF
    R2=csf(R,Z,viewingDis);
    A2=csf(A,Z,viewingDis);
    O2=csf(O,Z,viewingDis);

    % Saptial masking
    [R3,A3]=spatial_masking(R2,A2,Z,1);
    O3=abs(O2);
    
    % Temporal masking
    L3=O3-R3;
    if frIdx==1
        prevDWT=O2;
    end
    currDWT=O2;
    masker=abs(prevDWT-currDWT);
    L4=temporal_masking(L3,masker,0.5,Z);
    A4=temporal_masking(A3,masker,0.5,Z);
    O4=O3;
    
    % DLM and AIM
    framescoresDLM(frIdx)=Lpooling(L4,O4,Z,2,1,0.1); % 2: spatial pooling exponent 1: frequency pooling exponent 0.1: border extraction
    framescoresAIM(frIdx)=Apooling(A4,Z,2,1,0.1);

    % Prepare for the next frame
    prevDWT = O2;
    refFrameY_xn_minus_1=refFrameY_xn;
    disFrameY_xn_minus_1=disFrameY_xn;
    refFrameY_yn_minus_1=refFrameY_yn;
    disFrameY_yn_minus_1=disFrameY_yn;
end

% Temporal pooling
framescores=framescoresDLM+27.45*framescoresAIM;
framescoresProcessed(1)=framescores(1);
for i=2:size(framescores,2)
    delta=framescores(i)-framescoresProcessed(i-1);
    if delta<=0; 
        alpha=0.04; 
    else
        alpha=0.5; 
    end
    framescoresProcessed(i)=framescoresProcessed(i-1)+alpha*delta;
end
score=mean(framescoresProcessed);

end

function X2=csf(X,Z,viewingDis)
X2=X;
theta=pi*viewingDis*Z(end,1)/180;
totalLevel = size(Z,1)-2;
for level=1:totalLevel
    thetatmp = theta/2^(5-level);
    for ori=1:3
        c=dalycsf(thetatmp,ori);
        [st ed]=get_scale_Index(Z,level,ori);
        X2(st:ed)=c.*X(st:ed);
    end
end
end

function DLM=Lpooling(L,O,Z,spaPool,frPool,border)
L=abs(L);
O=abs(O);

%spatial pooling
OhArr = [];
LhArr = [];
totalLevel=size(Z,1)-2;
for level=1:totalLevel
    for orient=1:3
    [st,ed] = get_scale_Index(Z,level,orient);
    Oh = O(st:ed);
    Lh = L(st:ed);
    Oh = center(Oh,Z,level,border);
    Lh = center(Lh,Z,level,border);
    OhArr = [OhArr (sum(Oh.^spaPool))^(1/spaPool)];
    LhArr = [LhArr (sum(Lh.^spaPool))^(1/spaPool)];
    end
end

%frequency pooling
OhTotal = (sum(OhArr.^frPool))^(1/frPool);
LhTotal = (sum(LhArr.^frPool))^(1/frPool);
DLM = LhTotal/OhTotal;
end

function AIM=Apooling(A,Z,spaPool,frPool,border)
A=abs(A);

%spatial pooling
AhArr = [];
totalLevel=size(Z,1)-2;
for level=1:totalLevel
    for orient=1:3
    [st,ed] = get_scale_Index(Z,level,orient);
    Ah = A(st:ed);
    Ah = center(Ah,Z,level,border);
    AhArr = [AhArr (sum(Ah.^spaPool))^(1/spaPool)];
    end
end

%frequency pooling
AIM = (sum(AhArr.^frPool))^(1/frPool)/(numel(A));
end

function c=dalycsf(w,ori)
lamda=0.228;
Q = [pi/2 pi pi*0.75];
fQ=w/(0.15*cos(4*Q(ori))+0.85);
if w<3.4
    c=0.981;
else
    c=2.6*(0.0192+lamda*fQ)*exp(-(lamda*fQ)^1.1);
end
end

function [R2,A2]=spatial_masking(R,A,Z,m)
MR=cal_mask(R,Z,m);
MA=cal_mask(A,Z,m);
R2 = sub_mask(R,MA);
A2 = sub_mask(A,MR);
end

function MX=cal_mask(X,Z,m)
MX = zeros(size(X));
totalLevel = size(Z,1)-2;
for level=1:totalLevel
    [st1,ed1] = get_scale_Index(Z,level,1);
    w1 = X(st1:ed1);
    w1 = abs(w1);
    w1 = wfilter(w1,Z,level);
    [st2,ed2] = get_scale_Index(Z,level,2);
    w2 = X(st2:ed2);
    w2 = abs(w2);
    w2 = wfilter(w2,Z,level);
    [st3,ed3] = get_scale_Index(Z,level,3);
    w3 = X(st3:ed3);
    w3 = abs(w3);
    w3 = wfilter(w3,Z,level);
    
    w = m*(w1+w2+w3)./3;
    MX(st1:ed1) = w;
    MX(st2:ed2) = w;
    MX(st3:ed3) = w;
end

end

function X2 = sub_mask(X,M)
absX = abs(X);
X2 = max((absX - M),0);
end

function Maskee2=temporal_masking(Maskee,Masker,tm,Z)
Masker = mask(abs(Masker),Z);
Maskee=abs(Maskee);
Maskee2=max((Maskee-tm.*Masker),0);
end

function W=mask(W,S)
levelNo = size(S,1)-2;
for level=1:levelNo
    for i=1:3
    [st,ed] = get_scale_Index(S,level,i);
    w = W(st:ed);
    w = abs(w);
    w = wfilter(w,S,level);
    W(st:ed) = w;
    end
end
end

function  W = wfilter(W,S,level)
h=S(level+1,1);
w=S(level+1,2);
W = reshape(W,[h,w]);
W = filter2([1,1,1;1,2,1;1,1,1]/10,W);
W = W(:);
end

function W = center(W,S,level,edge)
h=S(level+1,1);
w=S(level+1,2);
W = reshape(W,[h,w]);
margH = round(max(1,h*edge));
margW = round(max(1,w*edge));
W = W(margH:h-margH+1,margW:w-margW+1);
W = W(:);
end
 
function [st,ed] = get_scale_Index(S,level,orient)
st = S(1,1)*S(1,2)+1;
ed = st+3*S(2,1)*S(2,2)-1;
for k=1:level-1
    st = ed + 1;
    ed = st+3*S(k+2,1)*S(k+2,2)-1;
end
    
coeffNo = (ed-st+1)/3;
st = st+(orient-1)*coeffNo;
ed = st+coeffNo-1;
end

function frameYUV = get_one_frame(seqAddr,frIdx,chrotype,height,width)
fid = fopen(seqAddr,'r');

frameYUV = {};
if strcmp(chrotype,'420')==1
    frameSize = height*width*(3/2);
    fseek(fid,frameSize*(frIdx-1),'bof');
    Y = fread(fid,[width height],'uint8=>single');
    U = fread(fid,[width/2 height/2],'uint8=>single');
    V = fread(fid,[width/2 height/2],'uint8=>single');
else
    display('error in get_one_frame().');
end

frameYUV{1} = Y';
frameYUV{2} = U';
frameYUV{3} = V';

fclose(fid);
end

function Y = crop_frame(Y,blocksize)
[height, width] = size(Y);
height2 = floor(height/blocksize)*blocksize;
width2 = floor(width/blocksize)*blocksize;
Y = Y(1:height2,1:width2);
end

function [R,A]=decouple4video(O,D,Z,totalLevel)
% approximation subband
st = 1;
ed = Z(1,1)*Z(1,2);
R(st:ed) = D(st:ed);

% other subbands
for level=1:totalLevel
    for i=1:3
    [st,ed] = get_scale_Index(Z,level,i);
    O1 = O(st:ed);
    D1 = D(st:ed);
    
    k=min(max(D1./(O1+1e-30),0),1);
    tmpH = k.*O1;
    
    R(st:ed) = tmpH;
    end
end

A = D-R;
end
