clear;
errorProbs = [1e-6 1e-5 1e-4 1e-3];
numSimulations = 5;
alamoutiErrors = zeros(length(errorProbs),1);
for i = 1:length(errorProbs)
    Pe = errorProbs(i);
    N = 100 * 1/Pe;
    errorRate = monteCarlo(@alamouti, N, numSimulations);
    alamoutiErrors(i) = errorRate;
end

figure
plot(errorProbs, alamoutiErrors, 'b--o');
title('Alamouti Code');
xlabel('Probability of Error');
ylabel('Symbol Error Rate');
saveas(gcf,'alamouti.png');

function [rate] = alamouti(N)
    snr = 10;
    M=4;%M-ary Constellation
    x = randi([0 1],N,1);
    y1 = qammod(x,M); % y1 is 1 x N
    y = zeros(2,N); % 2 x N 

    % beginning at column 1, set every other column to reshaped y1
    y(:,1:2:end) = reshape(y1,2,N/2);

    % beginning at column 2, set every other column to conjugate reshaped y1 and flip
    % y1* columns
    % element wise multiply y1*
    bitFlip = kron(ones(1,N/2),[-1;1]); % create 2 x N/2 matrix, one row is -1, one row is 1
    y(:,2:2:end) =(bitFlip.*flipud(reshape(conj(y1),2,N/2)));% [-x2* x1*; ....]

    h1 = 1/sqrt(2)*(randn(1,N/2) + 1i*randn(1,N/2)); % create 1 x N/2 vector of random complex numbers normally distributed
    h2 = 1/sqrt(2)*(randn(1,N/2) + 1i*randn(1,N/2));
    h=[h1;h2]; % h is 2 x N/2
    %disp(size(h));
    H = kron(h,ones(1,2)); % duplicate each column, place duplicate next to orig, H is 2 x N

    hEq = zeros(2,N); % 2 x N
    % same transformation done to y, set every other column of hEq to h
    hEq(:,(1:2:end)) = h; % [h1 0 ... ; h2 0...]
    hEq(:,(2:2:end)) = kron(ones(1,N/2),[1;-1]).*flipud(h); % [h1 h2 ..; h2 -h1..]
    hEq(1,:) = conj(hEq(1,:)); %  [h1* h2* ... ; h2 -h1 .... ] % take conj of row 1
    hEqPower = sum(hEq.*conj(hEq),1); % sum each column, 1 x N vector


    gaussianNoise=1/sqrt(2)*(randn(1,N)+1i*randn(1,N));% noise vector 1 x N
    % element wise multiply H by modulated signal, sum each column, add noise, 1 x N
    y3=sum(H.*y,1)+10^(-(snr-10*log10(20))/20)*gaussianNoise; %snr is in dB, y3 is 1 x N

    yMod = kron(reshape(y3,2,N/2),ones(1,2)); % duplicate each column, place duplicate next to orig, yMod is 2 x N
    yMod(2,:) = conj(yMod(2,:)); % take conjugate of second row
    yHat = sum(hEq.*yMod,1)./hEqPower; % [h1*y1 + h2y2*, h2*y1 -h1y2*, ... ]
    z=qamdemod(yHat,M);
    z=reshape(z,N,1);
    [num rate]=symerr(x,z);
end

function mcSymErroRate = monteCarlo(func, numSymbols, N)
    errSum = 0;
    for i = 1:N
        symErrorRate = func(numSymbols);
        %symErrorRate = getErrorRate(origSymbols, outSymbols);
        errSum = errSum + symErrorRate;
        fprintf('Done simulation %i, error rate %d\n', i, symErrorRate);
    end
    mcSymErroRate = errSum / N;
end