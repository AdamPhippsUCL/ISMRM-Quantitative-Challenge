%% Image settings

% Define TE
TE = 125;

% ROI
ROInum = 3;

%% Bootstrapping
Nboot = 50;


%% Iterate over b values

bvec = [500,1000,2000];

% Initialise arrays for results
ADCfits = zeros(3,1);
T2fits = zeros(3,1);
sigma0fits = zeros(3,1);

ADCboots = zeros(3,Nboot);
T2boots = zeros(3,Nboot);
sigma0boots = zeros(3,Nboot);

for bindx = 1:3

    b = bvec(bindx);
    
    % Load ROI values
    load(['C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Quantitative Challenge\Outputs\ROI Values\b' num2str(b) '\ROI' num2str(ROInum) '.mat'])

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
    nbin = 400;
    binnums = [400,500,800,1000];

    binedges = linspace(binmin, binmax, nbin+1);
    bincenters = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;
    
    
    f = figure('visible','off');
    H = histogram(ROIvals, binedges);
    hold on;
    counts = H.Values;
    Ncounts = sum(counts);
    close(f)


    % Check there are enough bins
    indx = 1;
    while and( max(counts)>0.1*Ncounts, indx < length(binnums))
        indx=indx+1;
        nbin = binnums(indx)
        binedges = linspace(binmin, binmax, nbin+1);
        bincenters = (binedges(:,1:end-1) + binedges(:,2:end)) / 2 ;
        f = figure('visible','off');
        H = histogram(ROIvals, binedges);
        hold on;
        counts = H.Values;
        close(f)
    end

    %% Fitting 

    F = figure;
    H = histogram(ROIvals, binedges);
    hold on;

    % Initial guess
    fdguess = mean(ROIvals);
    T2guess = 500;
    sigma0guess = 0.01;
    
    beta0guess = [sigma0guess, T2guess, fdguess];
    
    % Apply fitting
    [coeffs, resnorm] = fitDistToHist(counts, bincenters, beta0guess = beta0guess, TE = TE, disttype = 'Ratio');
    
    sigma0 = coeffs(1);
    T2 = coeffs(2);
    fd = coeffs(3);
    ADC = -log(fd)/b;

    ADCfits(bindx) = ADC;
    T2fits(bindx) = T2;
    sigma0fits(bindx) = sigma0;
    
 
    % Create and display distribution
    dz = (binmax-binmin)/nbin;
    [dist, signals] = RatioDistRician(exp(-TE/T2), exp(-TE/T2)*fd, sigma0, zmin = binmin, zmax = binmax, dz = dz, ymin=binmin, ymax=binmax, dy=dz);
    
    plot(signals, sum(counts)*dist*dz);



    %% Bootstrapping

    for bootindx = 1:Nboot
    
        Nsample = Ncounts;
    
        % Take random sample
        sample = ROIvals(randsample(Ncounts, Nsample));
    
        f = figure('visible','off');
        H = histogram(sample, binedges);
        hold on;
        counts = H.Values;
        close(f);
    
        % Initial guess
        fdguess = mean(sample);
        T2guess = 200 + 600*rand();
        sigma0guess = 0.005 + 0.02*rand();
        
        beta0guess = [sigma0guess, T2guess, fdguess];
        
        % Apply fitting
        [coeffs, resnorm] = fitDistToHist(counts, bincenters, beta0guess = beta0guess, TE = TE, disttype = 'Ratio');
        
        sigma0boots(bindx, bootindx) = coeffs(1);
        T2boots(bindx, bootindx) = coeffs(2);
        ADCboots(bindx, bootindx) = -log(coeffs(3))/b;
      
        % Create and display distribution
        dz = (binmax-binmin)/nbin;
        [dist, signals] = RatioDistRician(exp(-TE/coeffs(2)), exp(-TE/coeffs(2))*fd, coeffs(1),  zmin = binmin, zmax = binmax, dz = dz, ymin=binmin, ymax=binmax, dy=dz);
        plot(signals, Ncounts*dist*dz)
    
    
        pause(2)
    
    end


end

% Meta data
Meta = struct();
Meta.binmin = binmin;
Meta.binmax = binmax;
Meta.nbin = nbin;
Meta.Nboot = Nboot;
Meta.Ncounts = Ncounts;
Meta.Nsample = Nsample;


% Save results

OutputFolder = "C:\Users\adam\OneDrive - University College London\UCL PhD\PhD Year 1\Projects\ISMRM Quantitative Challenge\Outputs\Measurements";
OutputPath = [char(OutputFolder) '/ROI' num2str(ROInum)];

if ~exist(OutputPath, "dir")
   mkdir(OutputPath)
end

save([OutputPath '/ADCfits.mat'], 'ADCfits');
save([OutputPath '/T2fits.mat'], 'T2fits');
save([OutputPath '/sigma0fits.mat'], 'sigma0fits');

save([OutputPath '/ADCboots.mat'], 'ADCboots');
save([OutputPath '/T2boots.mat'], 'T2boots');
save([OutputPath '/sigma0boots.mat'], 'sigma0boots');

save([OutputPath '/Meta.mat'], 'Meta');