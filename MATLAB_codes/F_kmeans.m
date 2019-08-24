function [cluster_index,varargout] = F_kmeans(input_matrix,varargin)
%% [cluster_index,varargout] = F_kmeans(input_matrix,varargin)
%  F_kmeans performs unsuperized kmeans clustering with an optimal k. 
%  The number k is determined by Elbow methods. For more about Elbow method:
%  http://www.sthda.com/english/articles/29-cluster-validation-essentials/96-determining-the-optimal-number-of-clusters-3-must-know-methods/
%  You can also define the k for kmeans clustering. If doing so, the
%  funcion skips the step of calculating optimal k.
%
%  Example:
%  [cluster_IDX,PCA_score,~,~,tSNE_score] = F_kmeans(input_matrix,...
%                                                    'k_for_kmeans',5,...
%                                                    'SelectVariableGene','no',...
%                                                    'number_of_pc',50);
%  [cluster_IDX,PCA_score,PCA_coeff,PCA_explained,~] = F_kmeans(input_raw_matrix,...
%                                      'knnsmooth','yes',...
%                                      'k_for_knnsmooth',7,...
%                                      'PC_for_knnsmooth',15,...
%                                      'FTtransPCA','yes',...
%                                      'SelectVariableGene','no',...
%                                      'number_of_pc',50,...
%                                      'plot_PCA','yes',...
%                                      'plot_tSNE','no'...
%                                      );    
%
%   Pre-requirements: (should be in the same folder)
%       1) kmeans_opt.m (copyright@Sebastien De Landtsheer)
%               (https://www.mathworks.com/matlabcentral/fileexchange/65823-kmeans_opt-optimal-k-means-with-elbow-method)
%       2) knn_smooth.m (copyright@Florian Wagner, Maayan Baron)
%                    knn_smooth requires rsvd.m (randomized Singular Value Decomposition function)
%       3) cbrewer.m associated stuff in a folder(copyright@Charles Robert, http://colorbrewer.org/)
%               (https://www.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab)
%       4) Tested with Matlab R2017a based 'pca' and 'tsne' function.
%       5) F_normalize.m (copyright@Bo Xia)
%
%   input_matrix: gene_by_cell expression matrix. Type: double. Usually it
%                 is not normalized, but can be normalizaed matrix.
%   varargin:
%         'k_for_kmeans'        : k for doing kmeans clustering with self
%                                 assigned number of k. By doing so, the
%                                 function skips the the step of
%                                 calculating optimal k.
%                       Default : 0
%                       Optional: integers.
%         'knnsmooth'           : Do not do knn-smoothing.
%                       Default : 'no'
%                       Optional: 'yes'
%         'k_for_knnsmooth'     : k for doing knn-smoothing.
%                       Default : 3
%                       Optional: suggested as 2^(2:7)-1
%         'PC_for_knnsmooth'    : number of PCs for doing knn-smoothing.
%                       Default : 10
%                       Optional: suggest 10-50
%         'FTtransPCA'          : This is usually done after knn_smoothing.
%                                 The input matrix is normalized to median
%                                 and Freemen-Turkey transformed (which
%                                 keeps the variance of genes.FT-transform:
%                                           y = sqrt(x)+sqrt(x+1)
%                                 Then PCA will be done with all genes.
%                                 When 'FT_transPCA' is 'yes', default of
%                                 'SelectVariableGene' will be 'no'.
%                       Default : 'no'
%                       Optional: 'yes'
%         'SelectVariableGene'  : Do not select variable genes.
%                       Default : 'yes',fano_factor>1.5*nanmean(fano_factor).
%                                 If use 'yes', it assumes the input matrix
%                                 is raw expression matrix (but
%                                 pre-filtered bad-quality cells).                      
%                                 When 'FT_transPCA' is 'yes',
%                                 default of 'SelectVariableGene' will be
%                                 'no', unless specifically assign 'yes'.
%                       Optional: 'no'. In this case, the input matrix will
%                                 be used directly for PCA.
%         'number_of_pc'        : Number of PCs used in kmeans clustering and tSNE.
%                       Default : 50
%                       Optional: suggest integer of 10-100, can be 2-1000. 
%         'plot_PCA'            : plot 2D-scale embedding with PC1 and PC2, and
%                                 cumulative variance explained.
%                       Default : 'yes'
%                       Optional: 'no'
%         'plot_tSNE'           : plot 2D-scale embedding with tSNE1 and tSNE2.
%                                 Meanwhile, tSNE plot also contains
%                                 cluster annotations. Recommended!
%                       Default : 'yes'
%                       Optional: 'no'
%         'visual_mode'         : plot 2D-scale embedding with tSNE of PCA
%                                 with a bit random assignmend so that to
%                                 avoid overlapping embedding.
%                       Default : 'off'
%                       Optional: 'on'
%   varargout:
%     varargout{1} = 'PCA_score'; PCA score used for kmeans clustering and visualization;
%     varargout{2} = 'PCA_coeff'; PCA coefficients
%     varargout{3} = 'PCA_explained'; PCA variance explained
%     varargout{4} = 'tSNE_score'; tSNE score used for kmeans clustering and visualization
%
%   By Bo Xia, 2018-09-13

fprintf('************************************************************\n')
fprintf('Started kmeans clustering with parameters as below:\n')

%% Step 0: input parser

defaultSmooth = 'no';
validSmooth = {'no','yes'};
checkSmooth = @(x) any(validatestring(x,validSmooth));

defaultFTtransPCA = 'no';
validFTtransPCA = {'no','yes'};
checkFTtransPCA = @(x) any(validatestring(x,validFTtransPCA));

defaultSelectGene = 'yes';
validSelectGene = {'no','yes'};
checkSelectGene = @(x) any(validatestring(x,validSelectGene));

defaultPCAplot = 'yes';
validPCAplot = {'no','yes'};
checkPCAplot = @(x) any(validatestring(x,validPCAplot));

defaulttSNEplot = 'yes';
validtSNEplot = {'no','yes'};
checktSNEplot = @(x) any(validatestring(x,validtSNEplot));

defaultvisual_mode = 'off';
validvisual_mode = {'on','off'};
checkvisual_mode = @(x) any(validatestring(x,validvisual_mode));

defaultk  = 0;
defaultPC = 50;
defaultKnnSmoothK = 3;
defaultKnnSmoothPC = 10;

p = inputParser;
addOptional(p,'knnsmooth',defaultSmooth,checkSmooth)
addOptional(p,'FTtransPCA',defaultFTtransPCA,checkFTtransPCA)
addOptional(p,'SelectVariableGene',defaultSelectGene,checkSelectGene)
addOptional(p,'plot_PCA',defaultPCAplot,checkPCAplot)
addOptional(p,'plot_tSNE',defaulttSNEplot,checktSNEplot)
addOptional(p,'visual_mode',defaultvisual_mode,checkvisual_mode)
addParameter(p,'k_for_kmeans',defaultk,@isnumeric)
addParameter(p,'k_for_knnsmooth',defaultKnnSmoothK,@isnumeric)
addParameter(p,'PC_for_knnsmooth',defaultKnnSmoothPC,@isnumeric)
addParameter(p,'number_of_pc',defaultPC,@isnumeric)

parse(p,varargin{:})

defaultk = p.Results.k_for_kmeans;
if defaultk > 0
    fprintf(' Self-assigned k=%d for kmeans clustering. \n',defaultk)
else 
    fprintf(' Will calculate optimal k for kmeans clustering by Elbow method. \n')
end

defaultSmooth     = p.Results.knnsmooth;
disp([' knn-smoothing input matrix?  ',defaultSmooth])
if strcmp('yes',cellstr(defaultSmooth))
    defaultKnnSmoothK = p.Results.k_for_knnsmooth;
    disp(['   knnsmooth uses k =  ', num2str(defaultKnnSmoothK)])
    defaultKnnSmoothPC = p.Results.PC_for_knnsmooth;
    disp(['   knnsmooth uses PC_number = ', num2str(defaultKnnSmoothPC)])
end

defaultFTtransPCA = p.Results.FTtransPCA;
disp([' Use Freemen-Turkey transformation and then PCA?  ',defaultFTtransPCA])
if strcmp('yes',cellstr(defaultFTtransPCA))
    defaultSelectGene = 'no'; %re-assign this default when using FT transformation.
end

defaultSelectGene = p.Results.SelectVariableGene;
disp([' Select Variable Genes?  ',defaultSelectGene])

defaultPCAplot = p.Results.plot_PCA;
defaulttSNEplot = p.Results.plot_tSNE;
defaultvisual_mode = p.Results.visual_mode;

defaultPC         = p.Results.number_of_pc;
disp([' Number of PCs used for kmeans clustering:  ', num2str(defaultPC)])


%% Step 1: define input

%knnsmooth or not? dedault 'no'
if strcmp('no',cellstr(defaultSmooth))
    A = input_matrix;
elseif strcmp('yes',cellstr(defaultSmooth))
    A = knn_smooth(input_matrix,defaultKnnSmoothK,defaultKnnSmoothPC);
end 

%Freemen-Turkey transformation or not? dedault 'no'
if strcmp('yes',cellstr(defaultFTtransPCA))
    B = F_normalize(F_normalize(A,'median'),'FT');
elseif strcmp('no',cellstr(defaultFTtransPCA))
    B = A;
end 

%Select variable genes or not? default 'yes' or 'no' depends on condition
if strcmp('no',cellstr(defaultSelectGene))
    %dimention reduction input:
    d = B'; %transform the matrix to cell-by-gene, easier for PCA&tSNE
elseif strcmp('yes',cellstr(defaultSelectGene))
    %Select dynamic genes
    B_tpm = F_normalize(B,'tpm');
    B_log=log10(B_tpm+1);
    fano_factor  = var(B_tpm,0,2)./mean(B_tpm,2);
    %dimention reduction input:
    d = B_log(find(fano_factor>1.5*nanmean(fano_factor)),:)'; %transform
else
    disp '!!! Select Variable Genes?  Error!'
end
fprintf(' Number of genes used for PCA: %d \n',length(d(1,:)))
fprintf('********************\n')

%% Step 2: pca and kmeans clustering
[PCA_coeff, PCA_score, ~, ~, PCA_explained] = pca(d);

%kmeans clustering depends on if assigning number of k.
if defaultk == 0
    % k-means clustering with PCA scores.
    % This requires the kmeans_opt methods, which using ELBOW methods to
    % determine the optimum number of k.
    [cluster_index,centroids,SUMD,K] = kmeans_opt(PCA_score(:,1:defaultPC));
    fprintf(' Performed kmeans clustering with Elbow method-defined K = %d. \n',K)
elseif defaultk > 0
    K = defaultk;
    fprintf(' Performed kmeans clustering with assigned K = %d. \n',K)
    [cluster_index,centroids,SUMD] = kmeans(PCA_score(:,1:defaultPC),defaultk);
end

fprintf('Detected %d clusters; Number of each cluster is: \n',K)
for i=1: K
    fprintf('          Cluster %d:    %d \n',i,length(find(cluster_index==i)))
end

%% Step 3: visulization

%setup colormaps for clusting
cd ./cbrewer
if K >= 6
    type_cm = cbrewer('qual','Paired',K);
else
    type_cm = cbrewer('qual','Accent',5);
end
cd ..

if strcmp('on',cellstr(defaultvisual_mode))
    %add 1% random number for visualization
rng(1)
PCA_score = PCA_score + rand(length(PCA_score(:,1)),length(PCA_score(1,:))) * (abs(max(max(PCA_score)))+abs(min(min(PCA_score)))) ./ 100;
end

%PCA plot
if strcmp('yes',cellstr(defaultPCAplot))
    figure;
    gscatter(PCA_score(:,1),PCA_score(:,2),cluster_index,type_cm,'.',6);
    hold on;
        %Calculate mean of each cluster
        for i=1:K
            if length(find(cluster_index==i))>=2
            PCA_score_center(i,:) = mean(PCA_score(find(cluster_index==i),1:2));
            else
            PCA_score_center(i,:) = PCA_score(find(cluster_index==i),1:2);
            end
            text(PCA_score_center(i,1),PCA_score_center(i,2),num2str(i));
        end
    hold off;
    xlabel({'PC1, expalained variance: ',strcat(num2str(PCA_explained(1)),'%')});
    ylabel({'PC2, expalained variance: ',strcat(num2str(PCA_explained(2)),'%')});
    title (strcat('k-means (Elbow method) clustering and PCA embedding, k=',num2str(K)))
    
    figure;
    pc_max = length(PCA_score(1,:));
    if pc_max <= 100 & pc_max >= 30
    subplot(1,2,1);plot(cumsum(PCA_explained(1:20)),'-o');
    xlabel('PC 1-20'); ylabel('Cumulative PCA explained variance');
    subplot(1,2,2);plot(cumsum(PCA_explained(1:pc_max)),'-o');ylim([0 100])
    xlabel('PC 1-max'); 
    elseif pc_max > 100
    subplot(1,2,1);plot(cumsum(PCA_explained(1:20)),'-o');
    xlabel('PC 1-20'); ylabel('Cumulative PCA explained variance');
    subplot(1,2,2);plot(cumsum(PCA_explained(1:100)),'-o');ylim([0 100])
    xlabel('PC 1-100'); 
    elseif pc_max < 30
    plot(cumsum(PCA_explained(1:pc_max)),'-o');ylim([0 100])
    xlabel('PC 1-max'); ylabel('Cumulative PCA explained variance');
    end
    suptitle('Cumulative PCA explained variance');
    
end

%tSNEplot
    tSNE_score = tsne(d,'NumPCAComponents',defaultPC);
    
if strcmp('on',cellstr(defaultvisual_mode))
    %add 1% random number for visualization
rng(1)
tSNE_score = tSNE_score + rand(length(tSNE_score(:,1)),length(tSNE_score(1,:))) * (abs(max(max(tSNE_score)))+abs(min(min(tSNE_score)))) ./ 100;
end
    
if strcmp('yes',cellstr(defaulttSNEplot))
    figure;gscatter(tSNE_score(:,1),tSNE_score(:,2),cluster_index,type_cm,'.',6);
    hold on;
        %Calculate mean of each cluster
        for i=1:K
            if length(find(cluster_index==i))>=2
            tSNE_score_center(i,:) = mean(tSNE_score(find(cluster_index==i),:));
            else
            tSNE_score_center(i,:) = tSNE_score(find(cluster_index==i),1:2);
            end
            text(tSNE_score_center(i,1),tSNE_score_center(i,2),num2str(i));
        end
    hold off;
    xlabel('tSNE1');
    ylabel('tSNE2');
    title (strcat('k-means (Elbow method) clustering and tSNE embedding, k=',num2str(K)))

end

%% Step 4: outputs
    varargout{1} = PCA_score;
    varargout{2} = PCA_coeff;
    varargout{3} = PCA_explained;
    varargout{4} = tSNE_score;
%    varargout{5} = centroids;
    
fprintf('Finished kmeans clustering with %d clusters.\n',K)
if strcmp('yes',cellstr(defaultPCAplot)) || strcmp('yes',cellstr(defaulttSNEplot))
   disp('Check out the plot(s)!')
end
fprintf('************************************************************\n')
end
