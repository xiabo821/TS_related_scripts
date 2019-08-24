function [gene_order] = zavit_angle(M,angle)

M = zscore(M,0,2); 
[~, score] = pca(M); 
X = score(:,1); Y = score(:,2); 
zavit_angle = mod(atan2d(X,Y)+angle, 360) ;
[~, gene_order] = sort(zavit_angle); %gene order is the zavit ordering