clear

errorProbs = [1e-6 1e-5 1e-4 1e-3];
numSimulations = 5;
gaussianErrors = zeros(length(errorProbs),1);
rayErrors = zeros(length(errorProbs),1);
antSelErrors = zeros(length(errorProbs),1);
antSelErrors3 = zeros(length(errorProbs),1);
antSelErrors4 = zeros(length(errorProbs),1);
optComErrors = zeros(length(errorProbs),1);

for i = 1:length(errorProbs)
    Pe = errorProbs(i);
    numSymbols = 100 * 1/Pe;
    
    fprintf('Running MC for Gaussian\n');
    gaussSymErrorRate = monteCarlo(@gaussian, numSymbols, numSimulations, 0);
    gaussianErrors(i) = gaussSymErrorRate;

    fprintf('\nRunning MC for Rayleigh\n');
    raySymErrorRate = monteCarlo(@rayleigh, numSymbols, numSimulations, 0);
    rayErrors(i) = raySymErrorRate;

    fprintf('\nRunning MC for Antenna Selection\n');
    antSelSymErrorRate = monteCarlo(@antennaSelection, numSymbols, numSimulations, 2);
    antSelErrors(i) = antSelSymErrorRate;
    
    fprintf('\nRunning MC for Antenna Selection\n');
    antSelSymErrorRate = monteCarlo(@antennaSelection, numSymbols, numSimulations, 3);
    antSelErrors3(i) = antSelSymErrorRate;
    
    fprintf('\nRunning MC for Antenna Selection\n');
    antSelSymErrorRate = monteCarlo(@antennaSelection, numSymbols, numSimulations, 4);
    antSelErrors4(i) = antSelSymErrorRate;

    fprintf('\nRunning MC for Optimum Combining\n');
    optComSymErrorRate = monteCarlo(@optCombining, numSymbols, numSimulations, 0);    
    optComErrors(i) = optComSymErrorRate;
end

figure
plot(errorProbs, gaussianErrors, 'b--o');
title('Gaussian Channel');
xlabel('Probability of Error');
ylabel('Symbol Error Rate');
saveas(gcf,'gaussian.png');
clf

plot(errorProbs, rayErrors, 'b--o');
title('Rayleigh Channel');
xlabel('Probability of Error');
ylabel('Symbol Error Rate');
saveas(gcf,'rayleigh.png');
clf

plot(errorProbs, antSelErrors, 'b--o', 'DisplayName','2 antennas');
hold on
plot(errorProbs, antSelErrors3, 'r--o', 'DisplayName','3 antennas');
hold on
plot(errorProbs, antSelErrors4, 'g--o', 'DisplayName','4 antennas');
title('Antenna Selection');
xlabel('Probability of Error');
ylabel('Symbol Error Rate');
legend
saveas(gcf,'antsel.png');
clf

plot(errorProbs, optComErrors, 'b--o');
title('Optimum Combining');
xlabel('Probability of Error');
ylabel('Symbol Error Rate');
saveas(gcf,'optcom.png');
clf

function [rate] = gaussian(numSymbols, numAntennas)
    M = 4;
    %bitsPerSym = log2(M); % 4-QAM has 2 bits per symbol
    %data = randi([0 1],numSymbols*bitsPerSym,1); % numSymbols*bitsPerSym x 1 column vector of bits
    x = randi([0 1],numSymbols,1);    

    % bin is natural ordering of binary number
    tx = qammod(x,M,'bin','InputType','bit');
    %rx = awgn(tx,25); % adds white Gaussian noise to signal, second arg is dB
    awgnchan = comm.AWGNChannel;
    rx = awgnchan(tx);
    y = qamdemod(rx,M,'bin','OutputType','bit');
    
    [num rate]=symerr(x,y);
    %dataDec = getDecArr(data, bitsPerSym);
    %yDec = getDecArr(y, bitsPerSym);
end

function [rate] = rayleigh(numSymbols, numAntennas)
    M = 4;
    %bitsPerSym = log2(M); % 4-QAM has 2 bits per symbol
    x = randi([0 1],numSymbols,1); % numSymbols*bitsPerSym x 1 column vector of bits
    %symNoise = rand(numSymbols*bitsPerSym,1); % random numbers between 0 and 1
    %symNoise = randi([-1,1],numSymbols*bitsPerSym,1);
    %x = data + symNoise;
    %x(x>1)=1;
    %x(x<0)=0;

    % bin is natural ordering of binary number
    tx = qammod(x,M,'bin','InputType','bit');
    
    rayleighchan = comm.RayleighChannel();
    chanOut1 = rayleighchan(tx);
    
    y = qamdemod(chanOut1,M,'bin','OutputType','bit');
    [num rate]=symerr(x,y);
    %dataDec = getDecArr(data, bitsPerSym);
    %yDec = getDecArr(y, bitsPerSym);
end

function [rate] = antennaSelection(numSymbols, numAntennas)
    M = 4;
    %bitsPerSym = log2(M); % 4-QAM has 2 bits per symbol
    x = randi([0 1],numSymbols,1); % numSymbols*bitsPerSym x 1 column vector of bits

    % bin is natural ordering of binary number
    tx = qammod(x,M);
    channels = zeros(numSymbols, numAntennas);
    noiseArr = zeros(numSymbols, numAntennas);
    norms = zeros(numAntennas, 1);
    maxNormIdx = -1;
    maxNorm = -1;
    for i = 1:numAntennas
        h = 1/sqrt(2)*(randn(1,numSymbols) + 1i*randn(1,numSymbols));
        h = h.';
        channels(:,i) = h;
        tempNorm = norm(h);
        if tempNorm > maxNorm
            maxNormIdx = i;
            maxNorm = tempNorm;
        end
        norms(i) = tempNorm;
        noise = 1/sqrt(2)*(randn(1,numSymbols)+1i*randn(1,numSymbols));
        noise = noise.';
        noiseArr(:,i) = noise;
    end
    rx = tx.*(norms(maxNormIdx).*channels(:,maxNormIdx))+norms(maxNormIdx).*noiseArr(:,maxNormIdx);
    
    y1 = qamdemod(rx,M);
    [num1 rate]=symerr(x,y1);
    %dataDec = getDecArr(data, bitsPerSym);
    %yDec = getDecArr(y, bitsPerSym);
end


function [rate] = optCombining(numSymbols, numAntennas)
    M = 4;
    %bitsPerSym = log2(M); % 4-QAM has 2 bits per symbol
    x = randi([0 1],numSymbols,1); % numSymbols*bitsPerSym x 1 column vector of bits

    % bin is natural ordering of binary number
    tx = qammod(x,M);
    
    h1 = 1/sqrt(2)*(randn(1,numSymbols) + 1i*randn(1,numSymbols));
    h1 = h1.';
    h2 = 1/sqrt(2)*(randn(1,numSymbols) + 1i*randn(1,numSymbols));
    h2 = h2.';
    noise1=1/sqrt(2)*(randn(1,numSymbols)+1i*randn(1,numSymbols));
    noise1 = noise1.';
    noise2=1/sqrt(2)*(randn(1,numSymbols)+1i*randn(1,numSymbols));
    noise2 = noise2.';
    a1 = (abs(h1)).^2;
    a2 = (abs(h2)).^2;
    rx = tx.*(a1.*h1+a2.*h2)+a1.*noise1+a2.*noise2;
    
    y = qamdemod(rx,M);
    [num rate]=symerr(x,y);
    %dataDec = getDecArr(data, bitsPerSym);
    %yDec = getDecArr(y, bitsPerSym);
end

function decArr = getDecArr(binString, bitsPerSym)
    L = length(binString);
    tempBin = zeros(1, 2, 'uint8');
    decArrLen = L / bitsPerSym;
    decArr = zeros(1, decArrLen, 'uint8');
    decArrInd = 1;
    % matlab array indexes start at 1
    for i = 2:bitsPerSym:L
        tempBin(1) = binString(i-1);
        tempBin(2) = binString(i);
        dec = bi2de(tempBin);
        decArr(decArrInd) = dec;
        decArrInd = decArrInd + 1;
    end
end

function symErrorRate = getErrorRate(orig, output)
    diff = orig == output;
    diffsum = sum(diff(:) == 1);
    symErrorRate = diffsum / length(orig);
end

function mcSymErroRate = monteCarlo(func, numSymbols, N, numAntennas)
    errSum = 0;
    for i = 1:N
        symErrorRate = func(numSymbols, numAntennas);
        %symErrorRate = getErrorRate(origSymbols, outSymbols);
        errSum = errSum + symErrorRate;
        fprintf('Done simulation %i, error rate %d\n', i, symErrorRate);
    end
    mcSymErroRate = errSum / N;
end