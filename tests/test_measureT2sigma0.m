% MATLAB Script to measure T2 and sigma0 values from b=0 ROI values

%% Image settings

% Define TE
TE = 125;

% ROI
ROInum = 0;

% Load b=0 ROI values
b=0;
load(['C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Quantitative Challenge\Outputs\ROI Values\b' num2str(b) '\ROI' num2str(ROInum) '.mat'])

%% Fitting

% Define number of random permutes
Npermute = 5;

% Define histogram bins
binmin = 0.5;
binmax = 1.5;
nbin = 400;
binedges = linspace(binmin, binmax, nbin+1);
bincenters = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;

% Keep track of estimates
T2fits = zeros(Npermute, 1);
sigma0fits = zeros(Npermute, 1);

for pindx = 1:Npermute

    % Randomised permutation of values
    permuteROIvals = ROIvals(randperm(length(ROIvals)));

    % Distribution of ratios
    ratiovals = permuteROIvals./ROIvals; 

    % Display on histogram
    f= figure;%('visible','off');
    H = histogram(ratiovals, binedges);
    counts = H.Values;

    % Fit ratio distribution
    fdguess = 1;
    T2guess = 500;
    sigma0guess = 0.01;
    
    beta0guess = [sigma0guess, T2guess, fdguess];
    
    [coeffs, resnorm] = fitDistToHist(counts, bincenters, fd = 1, beta0guess = beta0guess, TE = TE, disttype = 'Ratio');

    sigma0 = coeffs(1)
    sigma0fits(pindx) = sigma0;

    T2 = coeffs(2)
    T2fits(pindx) = T2;

    dz = (binmax-binmin)/nbin;
    [dist, signals] = RatioDistRician(exp(-TE/T2), exp(-TE/T2), sigma0, zmin = binmin, zmax = binmax, dz = dz, ymin=binmin, ymax=binmax, dy=dz);
    hold on
    plot(signals, sum(counts)*dist*dz)
    % close(f);

end