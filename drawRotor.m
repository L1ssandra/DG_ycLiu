% drawRotor.m
Qrho = Q1;
Qu = Q2./Q1;
Qv = Q3./Q1;
QE = Q4;
QB1 = Q5;
QB2 = Q6;
gamma = 5/3;
QP = (gamma - 1)*(QE - 0.5*Qrho.*(Qu.^2 + Qv.^2) - 0.5*(QB1.^2 + QB2.^2));
QC = sqrt(abs(gamma*QP./Qrho));
QMach = sqrt(Qu.^2 + Qv.^2)./QC;
QBP = 0.5*(QB1.^2 + QB2.^2);

figure(1);
contourf(xc,yc,Q1,15);colormap(hot);
%mesh(xc,yc,Q1);colormap(cool);
title('Density')

figure(2);
contourf(xc,yc,QP,15);colormap(hot);
title('Pressure')

figure(3);
contourf(xc,yc,QBP,15);colormap(hot);
title('Magnetic pressure')

figure(4);
contourf(xc,yc,QMach,25);colormap(hot);
title('Mach number')
