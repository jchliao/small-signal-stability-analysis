% 读取 Matpower 数据
% mpc = loadcase('case39');
% disp(mpc);

function [U,I,Y] = pf(mpc)
% A = makeA(mpc.bus, mpc.branch);
myY = makeY(mpc.bus,mpc.branch);
[Y, Yf, Yt] = makeYbus(mpc.baseMVA, mpc.bus, mpc.branch);
S = getS(mpc.baseMVA,mpc.bus,mpc.gen);
Vsp = getVsp(mpc.bus,mpc.gen);
result = runpf(mpc);
printpf(result);

U = ones(size(Y,1),1);
for i = 1:10
    J = makeJ(Y,U,mpc.bus);
    b = makeDelta(S,Y,U,Vsp,mpc.bus);
    if norm(b, Inf) < 10e-6
        disp("收敛了")
        disp(i)
        Umag = abs(U);
        phase_deg = rad2deg(angle(U));
        U2 = [Umag,phase_deg];
        disp(U2)
        I = Y*U;
        break
    end
    dU = linsolve(J, b);
    dU = complex(dU(1:2:end, :),dU(2:2:end, :));
    U = U-dU;
end
end


function x = gauss(A, b)
[m, n] = size(A);
if m ~= n
    error('输入的矩阵不是方阵');
end
if numel(b) ~= m
    error('输入矩阵 A 和向量 b 的维度不匹配');
end
Ab = [A, b];
% 消元
for k = 1:m-1
    for i = k+1:m
        factor = Ab(i, k) / Ab(k, k);
        Ab(i, k:end) = Ab(i, k:end) - factor * Ab(k, k:end);
    end
end
% 回代
x = zeros(m, 1);
x(m) = Ab(m, end) / Ab(m, m);
for i = m-1:-1:1
    x(i) = (Ab(i, end) - Ab(i, i+1:m) * x(i+1:m)) / Ab(i, i);
end
end


function Delta = makeDelta(S,Y,U,Vsp,bus)
I = Y * U;
dS = S - U.*conj(I);
dP = real(dS);
dQ = imag(dS);
n = length(U);
Delta = zeros(2*n,1);
for i = 1:n
    type = bus(i,2);
    if type == 1
        Delta(2*i) = dP(i);
        Delta(2*i-1) = dQ(i);
    elseif type == 2
        dU2 = abs(Vsp(i)^2) - abs(U(i)^2);
        Delta(2*i) = dP(i);
        Delta(2*i-1) = dU2;
    else %type == 3
        dU = Vsp(i) - U(i);
        Delta(2*i-1) = real(dU);
        Delta(2*i) = imag(dU);
    end
end
end


function J = makesubJ(type,i,j,Y,I,U)
J = zeros(2,2);
if type == 3
    J(1,1) = -1;
    J(1,2) = 0;
    J(2,1) = 0;
    J(2,2) = -1;
    return;
end
G = real(Y(i,j));
B = imag(Y(i,j));
e = real(U(i));
f = imag(U(i));
a = real(I(i));
b = imag(I(i));
H = -(G*e+B*f);
N = B*e-G*f;
% J = N
L = -H;
if type == 1
    if i==j
        J(1,1) = H - a;
        J(1,2) = N - b;
        J(2,1) = N + b;
        J(2,2) = L - a;
    else
        J(1,1) = H;
        J(1,2) = N;
        J(2,1) = N;
        J(2,2) = L;
    end
elseif type == 2
    if i==j
        J(1,1) = H - a;
        J(1,2) = N - b;
        J(2,1) = -2*e;
        J(2,2) = -2*f;
    else
        J(1,1) = H;
        J(1,2) = N;
        J(2,1) = 0;
        J(2,2) = 0;
    end
end
J([1, 2], :) = J([2, 1], :);
end

%type=1是PQ节点，type=2是PV节点，type=3是平衡节点
function J = makeJ(Y,U,bus)
I = Y*U;
n = size(Y,1);
J = zeros(2*n,2*n);
for i = 1:n
    type = bus(i,2);
    for j = 1:n
        if type == 3
            J(2*i-1:2*i,2*i-1:2*i) = makesubJ(type,i,j,Y,I,U);
            break
        end
        if Y(i,j) == 0
            continue
        end
        J(2*i-1:2*i,2*j-1:2*j) = makesubJ(type,i,j,Y,I,U);
    end
end
end


% A 为 bus x line 的矩阵
function A = makeA(bus,line)
num_bus = size(bus,1);
num_line = size(line,1);
A = zeros(num_bus,num_line);
for i = 1:num_line
    fbus = line(i,1);
    tbus = line(i,2);
    A(fbus,i) = 1;
    A(tbus,i) = -1;
end
end

function Y = makeY(bus,line)
A = makeA(bus, line);
z = line(:,3) +  1j * line(:,4);
Y_0 = A * diag(1 ./ z) * A';
y_sh = 1 / 2 * 1j * line(:,5);
Y_sh = diag(diag(A * diag(y_sh) * A'));
Y = Y_0 + Y_sh;
end

function S = getS(baseMVA,bus,gen)
S = -(bus(:,3) + 1j*bus(:,4));
for i = 1:size(gen,1)
    Sg = gen(i,2) + 1j*gen(i,3);
    gen_bus = gen(i,1);
    S(gen_bus) = S(gen_bus) + Sg;
end
S = S/baseMVA;
end

function Vsp = getVsp(bus,gen)
[GEN_BUS, PG, QG, QMAX, QMIN, VG, MBASE, GEN_STATUS, PMAX, PMIN, ...
    MU_PMAX, MU_PMIN, MU_QMAX, MU_QMIN, PC1, PC2, QC1MIN, QC1MAX, ...
    QC2MIN, QC2MAX, RAMP_AGC, RAMP_10, RAMP_30, RAMP_Q, APF] = idx_gen;
Vm = bus(:,8);
Va = bus(:,9);
Va_rad = deg2rad(Va);
V = Vm .* exp(1i * Va_rad);
Vsp = V;
gen_bus = gen(:,GEN_BUS);
Vsp(gen_bus) = gen(:,VG);
end