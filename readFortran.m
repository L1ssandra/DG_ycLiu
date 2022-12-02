% drawtest.m

Xc = load('Xc.txt');
Yc = load('Yc.txt');

Q1 = load('Q1.txt');
Q2 = load('Q2.txt');
Q3 = load('Q3.txt');
Q4 = load('Q4.txt');
Q5 = load('Q5.txt');
Q6 = load('Q6.txt');

Nx = length(Xc);
Ny = length(Yc);

xc = zeros(Nx,Ny);
yc = zeros(Nx,Ny);

for j = 1:Ny
    xc(:,j) = Xc;
end

for i = 1:Nx
    yc(i,:) = Yc';
end

Q1 = reshape(Q1,Nx,Ny)';
Q2 = reshape(Q2,Nx,Ny)';
Q3 = reshape(Q3,Nx,Ny)';
Q4 = reshape(Q4,Nx,Ny)';
Q5 = reshape(Q5,Nx,Ny)';
Q6 = reshape(Q6,Nx,Ny)';

%drawRotor
%drawSmoothVortex
drawOTV