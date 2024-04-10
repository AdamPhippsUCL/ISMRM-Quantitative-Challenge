function [coeffs, resnorm] = fitDistToHist_combined(countsMatrix, bvec, bincenters, params, opts)

arguments

    countsMatrix % Matrix of counts (each row for each b value) [Nbval, Ncounts]
    bvec % vector of b values
    bincenters % array of bin centers

    % Distribution parameters
    params.sigma0
    params.T2
    params.ADC
    params.TE = 125
    params.Nav_ratio = [1,1,1]

    opts.disttype = 'Ratio'
    opts.beta0guess = [0.05, 500, 1e-3]

    opts.boundstype = 'fixed';
    opts.lb = [0.001, 250, 0.05e-3]
    opts.ub = [0.05, 1500, 2e-3]
    opts.boundsfrac = 0.25

end

%% STEPS

% 1. Normalise histogram counts
% 2. Define unknown parameter vector for fitting
% 3. Apply fitting


%% 1. Normalise histogram counts

RowSums = sum(countsMatrix, 2);
countsMatrix = countsMatrix./RowSums;


%% Define unknown parameters
Nparam = 3;


%% Fit


% Assign variable to base (for access in functions)
assignin('base', 'current_params', params)
assignin('base', 'bvec', bvec)
assignin('base', 'bincenters', bincenters)

% Define bounds
switch opts.boundstype

    case 'fixed'
        % Bounds set by user
        disp('');

    case 'dependent'
        opts.lb = (1-opts.boundsfrac)*opts.beta0guess;
        opts.ub = (1-opts.boundsfrac)*opts.beta0guess;


end



switch opts.disttype

    case 'Ratio'
        if Nparam == 3
            [coeffs, resnorm] = lsqcurvefit(@RatioDist3, opts.beta0guess(1:3), bincenters, countsMatrix, opts.lb(1:3), opts.ub(1:3) );
        end

end


end




function binfreqsMatrix = RatioDist3(b, x)

    % Generate ratio distribution given x
    % Outputs bin frequencies at bin centers predicted
    % for ratio distribution with sigma0
    
    arguments
        b % [sigma0, T2] value
        x % bin centers
    end
    
    % Get bin centers
    bincenters = evalin('base', 'bincenters');

    % Get bvec
    bvec = evalin('base', 'bvec');
    nbval = length(bvec);

    % Get parameters
    params = evalin('base', 'current_params');
    TE = params.TE;
    Nav_ratio = params.Nav_ratio;

    % Distribution parameters
    sigma0 = b(1);
    T2 = b(2);
    ADC = b(3);
    
    % b0 signal
    b0signal = exp(-TE/T2);
    
    % b signals (for each b value)
    bsignals = exp(-ADC*bvec)*b0signal;
    
    % Bin spacings
    nbin = length(bincenters);
    binmin = min(x);
    binmax = max(x);
    binspacing = x(2)-x(1);
    
    % Generate pdf over bin centers for each b value
    binfreqsMatrix = zeros(nbval, nbin);

    for bindx = 1:nbval

        pdfvals = RatioDistRician( ...
            b0signal, ...
            bsignals(bindx), ...
            sigma0, ...
            Nav_ratio=Nav_ratio(bindx), ...
            zmin = binmin, ...
            zmax = binmax, ...
            dz = binspacing,...
            ymin = binmin,...
            ymax = binmax,...
            dy = binspacing);
            
        % Evaluate bin frequencies
        binfreqs = pdfvals*binspacing;
        
        binfreqs = reshape(binfreqs, size(x));

        binfreqsMatrix(bindx,:) = binfreqs;

    end

end
