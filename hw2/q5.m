clear
y1 = zeros(2,1);
y1(1,1) = .5648 + 1i*1.6805;
y1(2,1) = .3930 - 1i*0.0210;

y2 = zeros(2,1);
y2(1,1) = .1632 + 1i*0.7613;
y2(2,1) = .2574 + 1i*1.1152;

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

symVecs = zeros(2,8);
symVecs(:,1) = 1/sqrt(2)*[1 + 1i; 1 + 1i];
symVecs(:,2) = 1/sqrt(2)*[1 + 1i; 1 - 1i];
symVecs(:,3) = 1/sqrt(2)*[1 - 1i; 1 + 1i];
symVecs(:,4) = 1/sqrt(2)*[1 - 1i; 1 - 1i];
symVecs(:,5) = 1/sqrt(2)*[-1 + 1i; -1 + 1i];
symVecs(:,6) = 1/sqrt(2)*[-1 + 1i; 1 - 1i];
symVecs(:,7) = 1/sqrt(2)*[1 - 1i; -1 + 1i];
symVecs(:,8) = 1/sqrt(2)*[-1 - 1i; -1 - 1i];

y1ML1 = norm(y1 - H1*symVecs(:,1));
y1ML2 = norm(y1 - H1*symVecs(:,2));
y1ML3 = norm(y1 - H1*symVecs(:,3));
y1ML4 = norm(y1 - H1*symVecs(:,4));
y1ML5 = norm(y1 - H1*symVecs(:,5));
y1ML6 = norm(y1 - H1*symVecs(:,6));
y1ML7 = norm(y1 - H1*symVecs(:,7));
y1ML8 = norm(y1 - H1*symVecs(:,8));

y2ML1 = norm(y2 - H2*symVecs(:,1));
y2ML2 = norm(y2 - H2*symVecs(:,2));
y2ML3 = norm(y2 - H2*symVecs(:,3));
y2ML4 = norm(y2 - H2*symVecs(:,4));
y2ML5 = norm(y2 - H2*symVecs(:,5));
y2ML6 = norm(y2 - H2*symVecs(:,6));
y2ML7 = norm(y2 - H2*symVecs(:,7));
y2ML8 = norm(y2 - H2*symVecs(:,8));