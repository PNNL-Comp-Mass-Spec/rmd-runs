rmd-runs
========

% A Multivariate Outlier Discovery Method for the Identification of Statistically 
% Extreme LC-MS Analyses to Improve Quality Control Processing of Proteomics Experiments.

% RMD_RUNS.m
% FEBUARY 2011
%
% DISCLAIMER
    % This computer software was prepared by Battelle Memorial Institute, hereinafter the Contractor, 
      under Contract No. DE-AC05-76RL0 1830 with the Department of Energy (DOE).  All rights in the 
      computer software are reserved by DOE on behalf of the United States Government and the Contractor 
      as provided in the Contract.  You are authorized to use this computer software for Governmental 
      purposes but it is not to be released or distributed to the public.  NEITHER THE GOVERNMENT NOR 
      THE CONTRACTOR MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF 
      THIS SOFTWARE.  This notice including this sentence must appear on any copies of this computer 
      software.

% CONTACT
    % Melissa Matzke, melissa.matzke@pnl.gov
    % Bobbie-Jo Webb-Robertson, bj@pnl.gov

% DESCRIPTION
    % The RMD_RUNS function will log10 transform peptide peak intensity,
    % that is, peptide abundance data and determine if any LC-MS analyses 
    % (ie,runs) in a N x P peptide data set are statistical outliers.  The 
    % statistical analysis is based on summarizing each LC-MS run as a vector 
    % of q=5 summary statistics which describe the peptide abundance distribution
    % for a specific run; a N x q matrix is then analyzed using robust PCA 
    % to compute a robust estimate of the covariance matrix used in a the
    % calculation of a robut Mahalanobis distance.

% REQUIREMENTS for COMPLETE FUNCTIONALITY
    % Statistics Toolbox
        % Obtain from MathWorks, Inc.
    % PNNL Specific Functions
        % Cell2NumericVector.m
        % CreateSequentialNumbers.m
        % CorrelationMatrix.m
        % r_mad.m
    % ECON Specific Functions (included in files)
	% robpca.m (Croux and Ruiz-Gazen, 2005)
	% L1median.m (Croux and Ruiz-Gazen, 2005)

% INSTALLATION
    % Extract .m files from zip folder and either
    	% Add destination file to MatLab Path
    	% Run in MatLab from destination folder


% REFERENCES
    % Robust PCA  
    	% Croux, C. and Ruiz-Gazen, A (2005), "High breakdown estimatros for 
    	% principal componets: the Projection-pursuit approach revisited", 
	% J Multivariate Analysis,95,206-226.
	% http://www.econ.kuleuven.be/public/NDBAE06/programs/#pca

