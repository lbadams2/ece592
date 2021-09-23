% lecture 6 slides 3-13
% part a
clear
H = zeros(2,2);
H(1,1) = 1 + .5i;
H(1,2) = -1 - .25i;
H(2,1) = .2 + .5i;
H(2,2) = -.3 - .7i;
[U,S,V] = svd(H);
N_s = 1;
f = V(:,N_s); % optimum beamformer
w = U(:,N_s); % optimum combiner


% part b
N_s = 1;
F = zeros(2,4);
F(:,1) = [1; 0];
F(:,2) = [0; 1];
F(:,3) = [1/sqrt(2); 1/sqrt(2)*1i];
F(:,4) = [1/sqrt(2); -1/sqrt(2)*1i];
HprodF = H*F; % columns are product of H with each codeword
norm1 = norm(HprodF(:,1));
norm2 = norm(HprodF(:,2));
norm3 = norm(HprodF(:,3));
norm4 = norm(HprodF(:,4)); % this is highest SNR