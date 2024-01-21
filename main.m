clear
casename = 'case9';
mpc = loadcase(casename); 
disp(mpc);
% mpc.gendata(:,10) = 0
[lam,A] = eig_it(mpc)

function [lamba,A] = eig_it(mpc)
[U,I,Y] = pf(mpc);
S = U.*conj(I);
S(abs(S)<1e-6) = 0;

Vx = real(U);
Vy = imag(U);

SL = (mpc.bus(:,3)+1j*mpc.bus(:,4))/mpc.baseMVA;
SG = S + SL;
IL = conj(SL./U);
IG = IL + I;

gendata = mpc.gendata;
% gendata(:,10) = 100
gennum = size(gendata,1);
genbus = gendata(:,1);

loadindex = abs(SL) > 1e-6;
% loadindex(genbus) = false;
loadbuslist = find(loadindex);

for i = 1:size(mpc.gen,1)
    genbus = gendata(i,1);
    ws = 2*pi*50;
    D  = gendata(i,10);
    Tj = gendata(i,2);
    Xd = gendata(i,5);
    Igen = IG(genbus);
    E = U(genbus) + 1j*Xd*Igen;
    v = [ws,-D/Tj,-1/Tj,-1,Vx(genbus),Vy(genbus),imag(E),-Xd,-real(E),Xd];
    r = [1,2,2,3,3,3,4,4,5,5];
    c = [2,2,3,3,4,5,1,5,1,4];
    gen_danamic(i).A = sparse(r,c,v,5,5);
    v = [imag(Igen),real(Igen),1,1];
    r = [3,3,4,5];
    c = [1,2,2,1];
    gen_danamic(i).b = sparse(r,c,v,5,2);
end


for i = 1:length(loadbuslist)
    loadbus = loadbuslist(i);
    V0 = abs(U(loadbus));
    Vx0 = Vx(loadbus);
    Vy0 = Vy(loadbus);
    Ix0 = real(IL(loadbus));
    Iy0 = imag(IL(loadbus));
    P0 = real(SL(loadbus));
    Q0 = imag(SL(loadbus));
    v = [1,-2/V0*P0,1,-2/V0*Q0,-V0,1,-Vx0,-Vy0,1,-Vy0,Vx0];
    r = [1,1,2,2,3,4,4,4,5,5,5];
    c = [1,3,2,3,3,1,4,5,2,4,5];
    load_danamic(i).A = sparse(r,c,v,5,5);
    v = [Vy0,Vx0,-Iy0,-Ix0,-Ix0,Iy0];
    r = [3,3,4,4,5,5];
    c = [1,2,1,2,1,2];
    load_danamic(i).b = sparse(r,c,v,5,2);
end

APart = blkdiag(gen_danamic(:).A,load_danamic(:).A);

YY = makeYY(Y);

for i = 1:length(gendata(:,1))
    bus = gendata(i,1);
    genB{i} = sparse(5,size(YY,2));
    genB{i}(:,bus*2-1:bus*2) = gen_danamic(i).b;
end
for i = 1:length(loadbuslist)
    bus = loadbuslist(i);
    loadB{i} = sparse(5,size(YY,2));
    loadB{i}(:,bus*2-1:bus*2) = load_danamic(i).b;
end
BPart = vertcat(genB{:},loadB{:});

genc = [0 0 0 1 0;0 0 0 0 1];
loadc = [0 0 0 -1 0;0 0 0 0 -1];
for i = 1:length(gendata(:,1))
    bus = gendata(i,1);
    genC{i} = sparse(size(YY,1),5);
    genC{i}(bus*2-1:bus*2,:) = genc;
end
for i = 1:length(loadbuslist)
    bus = loadbuslist(i);
    loadC{i} = sparse(size(YY,1),5);
    loadC{i}(bus*2-1:bus*2,:) = loadc;
end
CPart = horzcat(genC{:},loadC{:});
DPart = -YY;
APart(abs(APart)<1e-6)=0;
ABCD = [APart,BPart;CPart,DPart];
ABCD = full(ABCD);

fxindx = [1:5:5*size(gendata,1), 2:5:5*size(gendata,1)];
logical_indx = false(1, size(ABCD,1));
logical_indx(fxindx) = true;
ABCDr = [ABCD(logical_indx,:);ABCD(~logical_indx,:)];
ABCDrc = [ABCDr(:,logical_indx),ABCDr(:,~logical_indx)];

Ap = ABCDrc(1:2*gennum, 1:2*gennum);
Bp = ABCDrc(1:2*gennum, 2*gennum+1:end);
Cp = ABCDrc(2*gennum+1:end, 1:2*gennum);
Dp = ABCDrc(2*gennum+1:end, 2*gennum+1:end);

A = Ap - Bp*inv(Dp)*Cp;
% [xR,lamba,xL] = eig(A);
lamba = eig(A);
end

function YY = makeYY(Y)
YY = sparse(2*size(Y,1),2*size(Y,2));
YY(1:2:end,1:2:end) = -imag(Y);
YY(1:2:end,2:2:end) = real(Y);
YY(2:2:end,1:2:end) = real(Y);
YY(2:2:end,2:2:end) = imag(Y);
end
