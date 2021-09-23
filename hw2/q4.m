% Lecture 5 starting at slide 18
clear
H1 = zeros(2,2);
H1(1,1) = exp(1i*pi/3);
H1(1,2) = exp(1i*pi/2);
H1(2,1) = exp(1i*pi/7);
H1(2,2) = exp(-1i*pi/4);

H2 = zeros(2,2);
H2(1,1) = exp(1i*pi/6);
H2(1,2) = .3*exp(1i*pi);
H2(2,1) = exp(1i*pi/5);
H2(2,2) = .1*exp(1i*pi/8);

y1 = zeros(2,1);
y1(1,1) = .5648 + 1i*1.6805;
y1(2,1) = .3930 - 1i*0.0210;

y2 = zeros(2,1);
y2(1,1) = .1632 + 1i*0.7613;
y2(2,1) = .2574 + 1i*1.1152;

symVecs = zeros(2,8);
symVecs(:,1) = 1/sqrt(2)*[1 + 1i; 1 + 1i];
symVecs(:,2) = 1/sqrt(2)*[1 + 1i; 1 - 1i];
symVecs(:,3) = 1/sqrt(2)*[1 - 1i; 1 + 1i];
symVecs(:,4) = 1/sqrt(2)*[1 - 1i; 1 - 1i];
symVecs(:,5) = 1/sqrt(2)*[-1 + 1i; -1 + 1i];
symVecs(:,6) = 1/sqrt(2)*[-1 + 1i; 1 - 1i];
symVecs(:,7) = 1/sqrt(2)*[1 - 1i; -1 + 1i];
symVecs(:,8) = 1/sqrt(2)*[-1 - 1i; -1 - 1i];

% part a
sHat1 = inv(conj(H1)*H1)*conj(H1);
sHat2 = inv(conj(H2)*H2)*conj(H2);

% part b
N_o = .01; % variance
denom1 = inv(H1*conj(H1));
denom1 = N_o * denom1(1,1);
snr1 = 1 / denom1;

denom2 = inv(H2*conj(H2));
denom2 = N_o * denom2(2,2);
snr2 = 1 / denom2;

% part c
distance1 = norm(sHat1*y1 - symVecs(:,1));
distance2 = norm(sHat1*y1 - symVecs(:,2));
distance3 = norm(sHat1*y1 - symVecs(:,3));
distance4 = norm(sHat1*y1 - symVecs(:,4));
distance5 = norm(sHat1*y1 - symVecs(:,5));
distance6 = norm(sHat1*y1 - symVecs(:,6));
distance7 = norm(sHat1*y1 - symVecs(:,7));
distance8 = norm(sHat1*y1 - symVecs(:,8));

y2distance1 = norm(sHat2*y2 - symVecs(:,1));
y2distance2 = norm(sHat2*y2 - symVecs(:,2));
y2distance3 = norm(sHat2*y2 - symVecs(:,3));
y2distance4 = norm(sHat2*y2 - symVecs(:,4));
y2distance5 = norm(sHat2*y2 - symVecs(:,5));
y2distance6 = norm(sHat2*y2 - symVecs(:,6));
y2distance7 = norm(sHat2*y2 - symVecs(:,7));
y2distance8 = norm(sHat2*y2 - symVecs(:,8));