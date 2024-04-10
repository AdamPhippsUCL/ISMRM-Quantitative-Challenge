%% Image settings

% Define TE
TE = 125;

% ROI
ROInum =1;


%% Bootstrapping
Nboot = 50;



%% Read data from different b values

bvec = [1000];
nb = length(bvec);

% Define bins
binmin = 0;
binmax = 1;
nbin = 800;

binedges = linspace(binmin, binmax, nbin+1);
bincenters = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;

countsMatrix = zeros(nb, nbin);
figs = [];

for bindx = 1:nb

    b = bvec(bindx);
    
    % Load ROI values
    load(['C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Quantitative Challenge\Outputs\ROI Values\b' num2str(b) '\ROI' num2str(ROInum) '.mat'])
    
    if bindx == 1
        Fig = figure;
    end
    Hist = histogram(ROIvals, binedges);
    hold on;

    % Remove anomalies from ROIvals
    upperprctle = prctile(ROIvals, 75);
    lowerprcntle = prctile(ROIvals, 25);
    med = median(ROIvals);
    IQR = upperprctle-lowerprcntle;
    
    ROIvals = ROIvals(ROIvals<med+3*IQR);
    ROIvals = ROIvals(ROIvals>med-3*IQR);


    f = figure('visible','off');
    H = histogram(ROIvals, binedges);
    hold on;
    counts = H.Values;
    Ncounts = sum(counts);
    close(f)

    countsMatrix(bindx, :) = counts;

end


%% Apply distribution fitting

% Initial guess
ADCguess = 1e-3;
T2guess = 1000;
sigma0guess = 0.02;

beta0guess = [sigma0guess, T2guess, ADCguess];

% Apply fitting
[coeffs, resnorm] = fitDistToHist_combined(countsMatrix, bvec, bincenters, beta0guess = beta0guess, TE = TE, disttype = 'Ratio');

% Results
sigma0 = coeffs(1);
T2 = coeffs(2);
ADC = coeffs(3);

% Make distributions and display
for bindx = 1:nb
    b = bvec(bindx);
    dz = (binmax-binmin)/nbin;
    [dist, signals] = RatioDistRician(exp(-TE/T2), exp(-TE/T2)*exp(-b*ADC), sigma0, zmin = binmin, zmax = binmax, dz = dz, ymin=binmin, ymax=binmax, dy=dz);
    
    plot(signals, length(ROIvals)*dist*dz);
end