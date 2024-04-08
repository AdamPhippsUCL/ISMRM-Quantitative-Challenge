%% Image settings

% Define TE
TE = 125;

% b value
b=1000;

% ROI
ROInum = 5;

load(['C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Quantitative Challenge\Outputs\ROI Values\b' num2str(b) '\ROI' num2str(ROInum) '.mat'])


%% Histogram

% Remove anomalies from ROIvals
upperprctle = prctile(ROIvals, 75);
lowerprcntle = prctile(ROIvals, 25);
med = median(ROIvals);
IQR = upperprctle-lowerprcntle;

ROIvals = ROIvals(ROIvals<med+3*IQR);
ROIvals = ROIvals(ROIvals>med-3*IQR);

% Define bins
binmin = 0;
binmax = 1;
nbin = 500;

binedges = linspace(binmin, binmax, nbin+1);
bincenters = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;

f = figure;%('visible','off');
H = histogram(ROIvals, binedges);
hold on;
counts = H.Values;


%% Fitting

fdguess = mean(ROIvals);
T2guess = 500;
sigma0guess = 0.01;

beta0guess = [sigma0guess, T2guess, fdguess];

[coeffs, resnorm] = fitDistToHist(counts, bincenters, beta0guess = beta0guess, TE = TE, disttype = 'Ratio');

sigma0 = coeffs(1)
T2 = coeffs(2)
fd = coeffs(3)
ADC = -log(fd)/b


% Create and display distribution
dz = (binmax-binmin)/nbin;
[dist, signals] = RatioDistRician(exp(-TE/T2), exp(-TE/T2)*fd, sigma0, zmin = 0, zmax = 1, dz = dz, ymin=0, ymax=1, dy=dz);

plot(signals, sum(counts)*dist*dz)
