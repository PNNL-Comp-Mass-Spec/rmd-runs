function [keep_runs,suspect_runs,LOG2RMD,RMD_PVALUES] = RMD_RUNS(DATA,CLASS,TECH_REPS,str1,str1val,str2,str2val,str3,str3val,str4,str4val,str5,str5val)


% DESCRIPTION
    % This function will log10 transform peptide peak intensity, that is,
    % peptide abundance data and determines if any LC-MS analyses (ie,runs)
    % in a N x P peptide data set are statistical outliers.  The statistical
    % analysis is based on summarizing each LC-MS run as a vector of q=5
    % summary statistics which describe the peptide abundance distribution
    % for a specific run; a n x q matrix is then analyzed using robust PCA 
    % to compute a robust estimate of the covariance matrix used in a the
    % calculation of a robut Mahalanobis distance.
    
    % Three diagnostic plots are given as default output,
    % (1) Box plots of the peptide abundance distribution are displayed, 
    %     The ordering is based on TECH_REPS and CLASS.
    %     The group color is based on CLASS.
    % (2) The log2 robust Mahalanobis distance values are plotted,
    %     Statistical outlier runs are labeled.  
    %     Statistical outlier runs are those that sit above the red line,
    %     which represents the chi-square critical value.
    %     Blue downward triangles represent technical replicates of
    %     biological samples; Red downward triangles represent all
    %     technical replicates of a biological sample; and, Gray triangles
    %     below the red line are not statistical outliers.
    % (3) The run-by-run correlation plot are displayed.
    %     Statistical outlier runs are labeled.
    
    % Five optional diagnostic plots are provided upon request,
    % (4) Box plots of the peptide abundance distribution given in the
    % order of LC analysis,
    %     The bold number identifies the acutal LC analysis number (ie.
    %     "TRUE RUN" number).
    %     The horizontal number underneath the LC analysis number
    %     represents the run order id (which is determined by the ordering
    %     of CLASS and TECH_REP).
    % (5) Empirical PDF
    %     Statistical outlier runs are idenitified as red.
    % (6) Biplot PC1 vs. PC2
    % (7) Biplot PC1 vs. PC3
    %     NOTE:  Both biplots are returned from 'rPCA_biplot'
    % (8) Box plots of the robust PCA scores.

% INPUT
    % DATA (REQUIRED)
        % DATA is to be of size N x P where N is the number of LC analyses 
        % (runs) and P is the number of peptides.
        % DATA should NOT be given as LOG TRANSFORMED.
        % DATA should be given as raw peak intensity (peptide abundance)
        % values.
        
    % CLASS (REQUIRED)
        % CLASS is the unique categories in the data.  For example if there
        % are two factors: BMI with 2 levels and EXPOSURE with 3 levels,
        % there would be 6 unique classes.  CLASS must be of size n x 1, 
        % where n is the number of runs.  It is assumed that CLASS is a 
        % numeric vector of integer values or a cell vector of strings.

    % TECH_REPS (REQUIRED)
        % TECH_REPS is a vector of integer values that describe the 
        % replication for each unique  biological sample.  For example,
        % if there are 5 samples each with  3 technical replicates under 
        % the assumption biological replicates are next to each other in 
        % the data matrix, 
            %TECH_REPS = [1,1,1,2,2,2,3,3,3,4,4,4,5,5,5];
        % However, this variable does not have to be sequential in this
        % manner, but it cannot skip numbers.

    % OPTIONAL Arguments must be given as 'name',value, (e.g.,
        % 'pvalue_threshold',0.0001)

    % pvalue_threshold (OPTIONAL)
        % This is the statistical threshold used to select outlier runs.
        % The p-value threshold must be between 0 and 1. 
        % The default is set to 0.0001
        
    % true_run_order (OPTIONAL)
        % Vector of length n, with the order specified 1:n(t),
        % e.g., if samples in order of the data matrix are 1,2,3,4,5 were 
        % run in order 3,1,4,5,2 on the Mass spectrometer then 
        % TRUE_RUN_ORDER = [3,1,4,5,2];
        
    % empirical_pdf (OPTIONAL)
        % Set to 1 if this plot is desired
        
    % rPCA_biplot (OPTIONAL)
        % Set to 1 if this plot is desired
        
    % rPCA_boxplot (OPTIONAL)
        % Set to 1 if this plot is desired
    
        
% OUTPUT
    % keep_runs
        % A vector of run identification numbers (row numbers) to be
        % retained in the data set for further statistical analyses.

    % suspect_runs
        % A vector of run identification numbers, if there statistical
        % outliers, that are to be removed from the data set.
    
    %LOG2RMD
        % A N x 1 vector of log2 robust Mahalanobis distance values for
        % each LC run in the data set. These values can be compared to a
        % log2(chi-square(q,alpha)) critical value.
        
    %RMD_PVALUES
        % A N x 1 vector of p-values based on the comparison of the roubst
        % Mahalanobis distance to a chi-square(q,alpha) critical value.
               
% USAGE EXAMPLES
    % keep_runs = RMD_RUNS(DATA,CLASS,TECH_REPS);
    % [keep_runs,suspect_runs,LOG2RMD,RMD_PVALUES] = RMD_RUNS(DATA,CLASS,TECH_REPS,'pvalue_threshold',0.001,'rPCA_biplot',1);
    % [keep_runs,suspect_runs] = RMD_RUNS(DATA,CLASS,TECH_REPS,'empirical_pdf',1,'true_run_order',TRUE_RUN_ORDER);
    % [keep_runs,suspect_runs,LOG2RMD,RMD_PVALUES] = RMD_RUNS(DATA,CLASS,TECH_REPS,'pvalue_threshold',0.001,'true_run_order',TRUE_RUN_ORDER,'empirical_pdf',1,'rPCA_boxplot',1,'rPCA_biplot',1);
    
% OTHER FUNCTIONS REQUIRED
    % Statistics Toolbox
        % Obtain from MathWorks, Inc.
        
    % PNNL Specific Functions
        % Cell2NumericVector.m
        % CreateSequentialNumbers.m
        % CorrelationMatrix.m
        % r_mad.m
        % robpca.m
        % L1median.m

% REFERENCES
    % Robust PCA  
    % Croux, C. and Ruiz-Gazen, A (2005), "High breakdown estimatros for 
    % principal componets: the Projection-pursuit approach revisited", 
    % J Multivariate Analysis,95,206-226.
    % http://www.econ.kuleuven.be/public/NDBAE06/programs/#pca
    
% CONTACT
    % Melissa Matzke, melissa.matzke@pnl.gov
    % Bobbie-Jo Webb-Robertson, bj@pnl.gov
    
% REVISION HISTORY
    % Melissa Matzke & Bobbie-Jo Webb-Robertson 1/26/11:  Created function
    

    if nargin < 3
        error('At minimum 3 input is required: 1) Data matrix, 2) Class vector and 3) Tech Replicate vector')
    end
    if nargin == 4 || nargin == 6 || nargin == 8 || nargin == 10 || nargin == 12        
        error(['At minimum 3 input is required: 1) Data matrix, 2) Class vector and 3) Tech Replicate vector','Optional Arguments must be input as pairs', 'Optional Arguments are pvalue_threshold, true_run_order, empirical_pdf, rPCA_biplot, rPCA_boxplot'])
    end
    if nargin > 13
        error('Too many inputs into the function')
    end

    % LOG10 Data
    x=log10(DATA);
    [N P]=size(x);  % number of LC analyses (runs) x number of peptides
    
    % Data check of CLASS
    [tmp1,tmp2] = size(CLASS);
    if min(tmp1,tmp2) > 1
        error('CLASS must be a one-dimensional vector with length equal to the first dimension of DATA')
    end
    if max(tmp1,tmp2) ~= N
        error('CLASS must be a one-dimensional vector with length equal to the first dimension of DATA: Check to make sure DATA is not transposed')
    end
    if iscell(CLASS) == 1
        CLASS = Cell2NumericVector(CLASS);
    elseif isnumeric(CLASS) == 0
        error('CLASS must be a one-dimensional vector with length equal to the first dimension of DATA of type numeric or cell')
    end
    if length(unique(CLASS)) ~= max(CLASS)
        CLASS = CreateSequentialNumbers(CLASS);
    end
    
    % Data check of TECH_REPS
    [tmp1,tmp2] = size(TECH_REPS);
    if min(tmp1,tmp2) > 1
        error('TECH_REPS must be a one-dimensional vector with length equal to the first dimension of DATA')
    end
    if max(tmp1,tmp2) ~= N
        error('TECH_REPS must be a one-dimensional vector with length equal to the first dimension of DATA: Check to make sure DATA is not transposed')
    end
    if isnumeric(TECH_REPS) == 0
        error('TECH_REPS must be a one-dimensional vector with length equal to the first dimension of DATA of type numeric')
    end
    if length(unique(TECH_REPS)) ~= max(TECH_REPS)
        TECH_REPS = CreateSequentialNumbers(TECH_REPS);
    end

    OPTIONAL_ARGUMENTS = {'pvalue_threshold','true_run_order','empirical_pdf','rPCA_biplot','rPCA_boxplot'};
    % Give default optional arguments and then overwrite if they are given
    pvalue_threshold = 0.0001;
    TRUE_RUN_ORDER = 1:N;
    true_run_order_plot = 0;
    empirical_pdf = 0;
    rPCA_biplot = 0;
    rPCA_boxplot = 0;
    
    if nargin > 4
        % determine what str1 is
        a = strmatch(str1,OPTIONAL_ARGUMENTS);
        if isempty(a)
            error('The first argument must be one of the following five options: pvalue_threshold, true_run_order, empirical_pdf, rPCA_biplot, rPCA_boxplot')
        end
        if a == 1
            pvalue_threshold = str1val;
        elseif a == 2
            TRUE_RUN_ORDER = str1val;
            true_run_order_plot = 1;
        elseif a == 3
            empirical_pdf = str1val;
        elseif a == 4
            rPCA_biplot = str1val;
        else
            rPCA_boxplot = str1val;
        end
    end
    
    if nargin > 6
        a = strmatch(str2,OPTIONAL_ARGUMENTS);
        if isempty(a)
            error('The second argument must be one of the following five options: pvalue_threshold, true_run_order, empirical_pdf, rPCA_biplot, rPCA_boxplot')
        end
        if a == 1
            pvalue_threshold = str2val;
        elseif a == 2
            TRUE_RUN_ORDER = str2val;
            true_run_order_plot = 1;
        elseif a == 3
            empirical_pdf = str2val;
        elseif a == 4
            rPCA_biplot = str2val;
        else
            rPCA_boxplot = str2val;
        end
    end
    
    if nargin > 8
        a = strmatch(str3,OPTIONAL_ARGUMENTS);
        if isempty(a)
            error('The third argument must be one of the following five options: pvalue_threshold, true_run_order, empirical_pdf, rPCA_biplot, rPCA_boxplot')
        end
        if a == 1
            pvalue_threshold = str3val;
        elseif a == 2
            TRUE_RUN_ORDER = str3val;
            true_run_order_plot = 1;
        elseif a == 3
            empirical_pdf = str3val;
        elseif a == 4
            rPCA_biplot = str3val;
        else
            rPCA_boxplot = str3val;
        end
    end
    
    if nargin > 10
        a = strmatch(str4,OPTIONAL_ARGUMENTS);
        if isempty(a)
            error('The fourt argument must be one of the following five options: pvalue_threshold, true_run_order, empirical_pdf, rPCA_biplot, rPCA_boxplot')
        end
        if a == 1
            pvalue_threshold = str4val;
        elseif a == 2
            TRUE_RUN_ORDER = str4val;
            true_run_order_plot = 1;
        elseif a == 3
            empirical_pdf = str4val;
        elseif a == 4
            rPCA_biplot = str4val;
        else
            rPCA_boxplot = str4val;
        end        
    end
    
    if nargin > 12
        a = strmatch(str5,OPTIONAL_ARGUMENTS);
        if isempty(a)
            error('The fifth argument must be one of the following five options: pvalue_threshold, true_run_order, empirical_pdf, rPCA_biplot, rPCA_boxplot')
        end
        if a == 1
            pvalue_threshold = str5val;
        elseif a == 2
            TRUE_RUN_ORDER = str5val;
            true_run_order_plot = 1;
        elseif a == 3
            empirical_pdf = str5val;
        elseif a == 4
            rPCA_biplot = str5val;
        else
            rPCA_boxplot = str5val;
        end        
    end
        
    % Check each of the 5 optional variables
    if max(size(pvalue_threshold)) > 1 || isnumeric(pvalue_threshold) == 0 || pvalue_threshold < 0 || pvalue_threshold > 1
        error('The pvalue_threshold must be a single value between 0 and 1')
    end
    % Data check of TECH_REPS
    [tmp1,tmp2] = size(TRUE_RUN_ORDER);
    if min(tmp1,tmp2) > 1
        error('TRUE_RUN_ORDER must be a one-dimensional vector with length equal to the first dimension of DATA')
    end
    if max(tmp1,tmp2) ~= N
        error('TRUE_RUN_ORDER must be a one-dimensional vector with length equal to the first dimension of DATA')
    end
    if isnumeric(TECH_REPS) == 0
        error('TRUE_RUN_ORDER must be a one-dimensional vector with length equal to the first dimension of DATA of type numeric')
    end
    if length(unique(TRUE_RUN_ORDER)) ~= N
        error('TRUE_RUN_ORDER must have sequential numbers from 1 to N, no numbers can be repeated')
    end
    if max(size(empirical_pdf)) > 1 || isnumeric(empirical_pdf) == 0
        error('The empirical_pdf must be a single value of 0 for no plot and 1 to display the graph')
    end
    if empirical_pdf ~= 0 && empirical_pdf ~= 1
        error('The empirical_pdf must be a single value of 0 for no plot and 1 to display the graph')
    end
    if max(size(rPCA_biplot)) > 1 || isnumeric(rPCA_biplot) == 0
        error('The rPCA_biplot must be a single value of 0 for no plot and 1 to display the graph')
    end
    if rPCA_biplot ~= 0 && rPCA_biplot ~= 1
        error('The rPCA_biplot must be a single value of 0 for no plot and 1 to display the graph')
    end
    if max(size(rPCA_boxplot)) > 1 || isnumeric(rPCA_boxplot) == 0
        error('The rPCA_boxplot must be a single value of 0 for no plot and 1 to display the graph')
    end
    if rPCA_boxplot ~= 0 && rPCA_boxplot ~= 1
        error('The rPCA_boxplot must be a single value of 0 for no plot and 1 to display the graph')
    end
% =================== END OF DATA INGEST AND DATA CHECKING

    % Box Plot of each run.
    figure('Name','Box Plots of Peptide Abundance Distribtuion by LC Analysis')
    boxplot(x','plotstyle','compact','colorgroup',CLASS)

    print -djpeg 'Boxplots'

    % Corr - correlation between biological reps.
    [FullCorrMatrix,Corr] = CorrelationMatrix(x,CLASS);

    % CALCULATE SKEWNESS for EACH RUN.
    Skewness_runs=skewness(x,0,2);

    % CALCULATE KURTOSIS for EACH RUN.
    Kurtosis_runs=kurtosis(x,0,2)-3;

    % CALCULATE FRACTION MISSING (NaN) VALUES FOR EACH RUN.
    f_missing_runs=sum(isnan(x'))./P;
    f_missing_runs=f_missing_runs(:);

    % CALCULATE MEDIAN ABSOLUTE DEVIATION (MAD) for the nth LC-MS RUN.
    mad_runs=mad(x,1,2);  

    % CREATE METRIC MATRIX
    % N runs x 5 quality measures
    Quality_Matrix = [f_missing_runs, mad_runs,...
                        Kurtosis_runs,Skewness_runs,Corr];
                
    % ROBUST PCA
    q=size(Quality_Matrix,2);
    [lambda,eigvecs,scores]=robpca(Quality_Matrix,q,@r_mad);

    % Percent Explained Variation
    P_Explained_Variation=lambda/sum(lambda);

    % Robust Covariance Estimate
    C_Sn=zeros(q,q);
    eigvecs_t=eigvecs';
    for jj=1:q
        C_Sn=C_Sn+(lambda(jj)*(eigvecs(:,jj)*eigvecs_t(jj,:)));
    end
    
    % ROBUST Mahalanobis Distance
    QM_Median=median(Quality_Matrix,1);
    dM_r=NaN(N,1);
    for i=1:N
        dM_r(i)= (Quality_Matrix(i,:)-QM_Median)*inv(C_Sn)*(Quality_Matrix(i,:)-QM_Median)';
    end
    LOG2RMD = log2(dM_r);
    RMD_PVALUES = 1 - chi2cdf(dM_r,q);
 
    suspect_runs = find(RMD_PVALUES <= pvalue_threshold);
    suspect_runs_dM=dM_r(suspect_runs)';
    keep_runs = find(RMD_PVALUES > pvalue_threshold);
    keep_runs_dM=dM_r(keep_runs)';
    
    % Match suspect runs to biological samples.
    Outlier_Runs = suspect_runs;
    Outlier_Samples_IDX = [];
    Outlier_Tech_rep_IDX = [];
 
    n = length(TECH_REPS);
    NS = max(TECH_REPS);
 
    a = zeros(n,1);
    a(Outlier_Runs) = 1;
    a = a.*TECH_REPS;
 
    for i = 1:NS
        count1 = length(find(TECH_REPS == i));
        count2 = length(find(a == i));
        if count2 > 0
            if count1 == count2
                % sample outlier
                Outlier_Samples_IDX = [Outlier_Samples_IDX;find(a == i)];
            else
                Outlier_Tech_rep_IDX = [Outlier_Tech_rep_IDX;find(a == i)];
            end
        end    
    end

    outlier_samples_dM=dM_r(Outlier_Samples_IDX)';
    outlier_tech_reps_dM=dM_r(Outlier_Tech_rep_IDX)';


% ============ Plot log2 Mahanalobis distance (dM_r) ============ %
figure('Name','Robust Mahanalobis Distance')
plot(keep_runs,log2(keep_runs_dM),'^','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor',[.5 .5 .5])
hold on
xlim([-2 N+3])
plot(Outlier_Samples_IDX,log2(outlier_samples_dM),'vr','MarkerFaceColor','r')
plot(Outlier_Tech_rep_IDX,log2(outlier_tech_reps_dM),'vb','MarkerFaceColor','b')
set(gca,'xtick',suspect_runs)
% Make new tick labels
xlabs=num2str(suspect_runs);
set(gca,'XTickLabel',{xlabs});
% text(suspect_runs,repmat(min(log2(dM_r))-1,length(suspect_runs),1),xlabs,...
%     'HorizontalAlignment','right','rotation',90,'FontSize',6)

pv_thres = 1 - pvalue_threshold;
% Draw horizontal red line for critical value
line([-2 N+3],[log2(chi2inv(pv_thres,q)) log2(chi2inv(pv_thres,q))],...
    'LineWidth',2,'Color','r')
% Label outlier runs
text(suspect_runs,log2(suspect_runs_dM)+0.1,xlabs,'FontSize',8)

ylabel('log2(Robust Mahalanobis Distance)','FontSize',10)
xlabel({'','LC-MS Identification Number'},'FontSize',10)
hold off;

print -djpeg 'rMdPlot'
% ============ Plot log2 Mahanalobis distance (dM_r) ============ %

% ============ Correlation Plot with Outlier Runs Labeled. ============ %
figure('Name','Correlation Heat Map');
imagesc(FullCorrMatrix)
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
colormap hot
colorbar('location','EastOutside')
caxis([0 1])
hold on;
text(suspect_runs,repmat(N+1,length(suspect_runs),1),xlabs,...
    'HorizontalAlignment','right','rotation',90,'FontSize',6)
text(zeros(length(suspect_runs),1),suspect_runs,xlabs,...
    'HorizontalAlignment','right','rotation',0,'FontSize',6)
hold off;

print -djpeg 'CorrelationHeatMap'
% ============ Correlation Plot with Outlier Runs Labeled. ============ %

    
    
%%%%%%%%%%%%%%%%%%OPTIONAL PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% ============ TRUE RUN ORDER BOX PLOT ============%
    if true_run_order_plot == 1
        [~, b]=sort(TRUE_RUN_ORDER);
        y=x(b,:);
        
        [c d]=sort(TRUE_RUN_ORDER(suspect_runs));
        suspect_true_runs=c;
        sr=suspect_runs(d);
      
        figure('Name','Box Plots Peptide Abundance by LC Run Order')
        boxplot(y','plotstyle','compact')
        set(gca,'XTickLabel',{' '})
        set(gca,'xtick',suspect_true_runs)
        % Label outlier LC runs  
        xlabs=num2str(suspect_true_runs);   % By LC Run Order
        text(suspect_true_runs,repmat(nanmin(nanmin(y))-.5,length(suspect_true_runs),1),xlabs,...
            'HorizontalAlignment','right','rotation',90,'FontSize',9,'FontWeight','bold')
        
        xlabs=num2str(sr);                  % By LC Run ID number
        text(suspect_true_runs,repmat(nanmin(nanmin(y))-1,length(suspect_true_runs),1),xlabs,...
            'HorizontalAlignment','center','rotation',0,'FontSize',8)

       print -djpeg 'Boxplot_TrueRunOrder'
    end
    % ============ TRUE RUN ORDER BOX PLOT ============%


% ============ Plot each run as Normal pdf. ============ %
    if empirical_pdf == 1
        % SUSPECT Runs in RED.
        figure('Name','Empirical PDF')
        for i=1:N
            x_i=sort(x(i,:));
            p_i=normpdf(x_i,nanmean(x_i),nanstd(x_i));
            plot(x_i,p_i,'LineStyle','-','Color',[0.5 0.5 0.5]); hold on;
            %ksdensity(x_i)
        end;
        s=x(suspect_runs,:);
        for i=1:length(suspect_runs)
            s_i=sort(s(i,:));
            p_i=normpdf(s_i,nanmean(s_i),nanstd(s_i));
            plot(s_i,p_i,'LineStyle','-','Color','r'); hold on;
        end;
        hold off;

        print -djpeg 'EmpiricalPDF'
    end
    % ============ Plot each run as Normal pdf. ============ %



% ============ Biplot rPCA results ============ %
    if rPCA_biplot == 1
        vbls = {'Fraction Missing','MAD','Kurtosis','Skew','R'};

    figure('Name','Biplot PC1 vs. PC2')
    biplot(eigvecs(:,[1,2]),'scores',scores(:,[1,2]),'varlabels',vbls);
    xlabel(['Component 1 (',num2str(P_Explained_Variation(1)*100,'%3.1f'),'%)'],'FontSize',11);
    ylabel(['Component 2 (',num2str(P_Explained_Variation(2)*100,'%3.1f'),'%)'],'FontSize',11)

    print -djpeg 'BiplotPC1PC2'

    figure('Name','Biplot PC1 vs. PC3')
    biplot(eigvecs(:,[1,3]),'scores',scores(:,[1,3]),'varlabels',vbls);
    xlabel(['Component 1 (',num2str(P_Explained_Variation(1)*100,'%3.1f'),'%)'],'FontSize',11)
    ylabel(['Component 3 (',num2str(P_Explained_Variation(3)*100,'%3.1f'),'%)'],'FontSize',11)

    print -djpeg 'BiplotPC1PC3'
    end
    % ============ Biplot rPCA results ============ %



% ============ Boxplot rPCA results ============ %
    if rPCA_boxplot == 1
        figure('Name','Boxplot of rPCA Scores by Principal Component')
        boxplot(scores)
        xlabel('Principal Component','FontSize',10)
        ylabel('Scores','FontSize',10)

        print -djpeg 'Boxplot_rPCAScores'
    end
    % ============ Boxplot rPCA results ============ %




