function [MAT]=F_normalize(mat,norm_method)
%%[MAT]=F_normalize(mat,norm_method)
%F_normalize requires input matrix [mat] and normlization methods(norm_method).
%
%[mat]: input gene expression matrix 
%   --genes as the rows, cells/samples as the columns
%   --usually UMI count matrix
%   
%(norm_method): available normlization methods:
%   'Anscombe', indicates Anscombe transformation
%   'mean', indicates mean-normalization
%   'median', indicates median-normalization
%   'tpm' --> normalization to 100000, but not 1000000
%   'tpm_log2' -->normalization to 100000, followed by log2 transformation
%   'tpm_log10' -->normalization to 100000, followed by log10 transformation
%   'FT', indicates Freeman_Turkey transformation
%
% Bo Xia

%%normalize

norm={'Anscombe','mean','median','tpm','tpm_log2','tpm_log10','FT'};

if string(norm_method)==norm(1); % Anscombe
    MAT=sqrt(mat+3/8)*2;
    fprintf('Done with %s transformation\n',string(norm(1)));

elseif string(norm_method)==norm(2); % mean
    MAT=mean(sum(mat))*bsxfun(@rdivide,mat,sum(mat));
    fprintf('Done with %s normalization\n',string(norm(2)));

elseif string(norm_method)==norm(3); % median
    MAT=median(sum(mat))*bsxfun(@rdivide,mat,sum(mat));
    fprintf('Done with %s normalization\n',string(norm(3)));

elseif string(norm_method)==norm(4); % tpm
    MAT=100000*bsxfun(@rdivide,mat,sum(mat));
    fprintf('Done with %s normalization\n',string(norm(4)));

elseif string(norm_method)==norm(5); % tpm_log2
    MAT=log2(1+100000*bsxfun(@rdivide,mat,sum(mat)));
    fprintf('Done with %s transformation\n',string(norm(5)));

elseif string(norm_method)==norm(6); % tpm_log10
    MAT=log10(1+100000*bsxfun(@rdivide,mat,sum(mat)));
    fprintf('Done with %s transformation\n',string(norm(6)));
    
elseif string(norm_method)==norm(7); % Freeman-Turkey
    MAT=sqrt(mat)+sqrt(mat+1);
    fprintf('Done with Freeman-Turkey transformation\n');
    
else
    fprintf('Your transformation, ‘%s’ method is not qualified.\n',norm_method);
    fprintf('Qualified methods listed as below: \n ')
        norm'
    
end