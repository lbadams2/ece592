clear
Pe = .001;
numSymbols = 100 * 1/Pe;
numSimulations = 20;
fprintf('Running MC for Gaussian\n');
gaussSymErrorRate = monteCarlo(@gaussian, numSymbols, 20);

fprintf('Running MC for Rayleigh\n');
raySymErrorRate = monteCarlo(@rayleigh, numSymbols, numSimulations);

function [dataDec, yDec] = gaussian(numSymbols)
    M = 4;
    bitsPerSym = log2(M); % 4-QAM has 2 bits per symbol
    data = randi([0 1],numSymbols*bitsPerSym,1); % numSymbols*bitsPerSym x 1 column vector of bits
    %symNoise = rand(numSymbols*bitsPerSym,1); % random numbers between 0 and 1
    symNoise = randi([-1,1],numSymbols*bitsPerSym,1);
    x = data + symNoise;
    x(x>1)=1;
    x(x<0)=0;

    % bin is natural ordering of binary number
    tx = qammod(x,M,'bin','InputType','bit');
    rx = awgn(tx,25); % adds white Gaussian noise to signal, second arg is dB
    y = qamdemod(rx,M,'bin','OutputType','bit');
    
    dataDec = getDecArr(data, bitsPerSym);
    yDec = getDecArr(y, bitsPerSym);
end

function [dataDec, yDec] = rayleigh(numSymbols)
    M = 4;
    bitsPerSym = log2(M); % 4-QAM has 2 bits per symbol
    data = randi([0 1],numSymbols*bitsPerSym,1); % numSymbols*bitsPerSym x 1 column vector of bits
    %symNoise = rand(numSymbols*bitsPerSym,1); % random numbers between 0 and 1
    symNoise = randi([-1,1],numSymbols*bitsPerSym,1);
    x = data + symNoise;
    x(x>1)=1;
    x(x<0)=0;

    % bin is natural ordering of binary number
    tx = qammod(x,M,'bin','InputType','bit');
    
    rayleighchan = comm.RayleighChannel(...
    'SampleRate',10e3, ...
    'PathDelays',[0 1.5e-4], ...
    'AveragePathGains',[2 3], ...
    'NormalizePathGains',true, ...
    'MaximumDopplerShift',30, ...
    'DopplerSpectrum',{doppler('Gaussian',0.6),doppler('Flat')}, ...
    'RandomStream','mt19937ar with seed', ...
    'Seed',22, ...
    'PathGainsOutputPort',true);
    [chanOut1,pathGains1] = rayleighchan(tx);
    
    y = qamdemod(chanOut1,M,'bin','OutputType','bit');
    
    dataDec = getDecArr(data, bitsPerSym);
    yDec = getDecArr(y, bitsPerSym);
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

function mcSymErroRate = monteCarlo(func, numSymbols, N)
    errSum = 0;
    for i = 1:N
        [origSymbols, outSymbols] = func(numSymbols);
        symErrorRate = getErrorRate(origSymbols, outSymbols);
        errSum = errSum + symErrorRate;
        fprintf('Done simulation %i, error rate %d\n', i, symErrorRate);
    end
    mcSymErroRate = errSum / N;
end