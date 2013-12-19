function [CM,CM_Class1,CM_Class2,CM_Class3,CM_Class4] = CorrelationMatrix(X_MATRIX,Class1,Class2,Class3,Class4)

%Author:  Bobbie-Jo Webb-Robertson
%Date Created:  05/25/2010
%Date Modified: 05/25/2010

% Revision History
    % Author:   Date:
    % revisions made

% This function takes a matrix (with missing data), and creates the
% correlation matrix

% Input
    % X:            An N x P matrix where N are the number of runs and P is
                    % the number of peptides

% Optional Input
    % Class1:       A (nx1) vector or cell structure that contains the class
                    % membership of each run to be averaged
    % Class2:       A (nx1) vector or cell structure that contains the class
                    % membership of each run to be averaged
    % Class3:       A (nx1) vector or cell structure that contains the class
                    % membership of each run to be averaged
    % Class4:       A (nx1) vector or cell structure that contains the class
                    % membership of each run to be averaged

% Output
    % CM:           the N x N correlation matrix
    % CM_Class      the average correlation within defined classes
                    % (ignoring self)
    % CM_TR:        the average correlation between technical replicates (ignoring self).
        
    
    
    
    
    
    if nargin < 1
        error('At least the matrix (runs x peptides) must be supplied')
    end
    
    if nargin > 1
        if iscell(Class1) == 1
            Class1 = Cell2NumericVector(Class1);
        elseif isnumeric(Class1) == 0
            error('Class Vectors must be either numeric or cell structures')
        end
    end

    if nargin > 2
        if iscell(Class2) == 1
            Class2 = Cell2NumericVector(Class2);
        elseif isnumeric(Class2) == 0
            error('Class Vectors must be either numeric or cell structures')
        end
    end

    if nargin > 3
        if iscell(Class3) == 1
            Class3 = Cell2NumericVector(Class3);
        elseif isnumeric(Class3) == 0
            error('Class Vectors must be either numeric or cell structures')
        end
    end

    if nargin > 4
        if iscell(Class4) == 1
            Class4 = Cell2NumericVector(Class4);
        elseif isnumeric(Class4) == 0
            error('Class Vectors must be either numeric or cell structures')
        end
    end

    [N,P] = size(X_MATRIX);
    CM = ones(N,N);
    for i = 1:N
        for j = i+1:N
            nancov(X_MATRIX(i,:),X_MATRIX(j,:));
            CM(i,j) = ans(1,2)/sqrt(ans(1,1)*ans(2,2));
            CM(j,i) = CM(i,j);
        end
    end
    
    if nargin > 1
        CM_Temp = CM;
        for i = 1:N
            CM_Temp(i,i) = NaN;
        end
        
        CM_Class1 = zeros(N,1);
        for i = 1:N
            a = find(Class1 == Class1(i));
            CM_Class1(i) = nanmean(CM_Temp(i,a));
        end
        z = isnan(CM_Class1);
        CM_Class1(z == 1) = 1;
    end
                    
    if nargin > 2
        CM_Class2 = zeros(N,1);
        for i = 1:N
            a = find(Class2 == Class2(i));
            CM_Class2(i) = nanmean(CM_Temp(i,a));
        end
        z = isnan(CM_Class2);
        CM_Class2(z == 1) = 1;
    end

    if nargin > 3
        CM_Class3 = zeros(N,1);
        for i = 1:N
            a = find(Class3 == Class3(i));
            CM_Class3(i) = nanmean(CM_Temp(i,a));
        end
        z = isnan(CM_Class3);
        CM_Class3(z == 1) = 1;
    end

    if nargin > 4
        CM_Class4 = zeros(N,1);
        for i = 1:N
            a = find(Class4 == Class4(i));
            CM_Class4(i) = nanmean(CM_Temp(i,a));
        end
        z = isnan(CM_Class4);
        CM_Class4(z == 1) = 1;
    end
