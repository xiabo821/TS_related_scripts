Q1_load_data=0;
colormaps=1;
Q2_conserve=0;
Q3_joint_PCA=0;
plot_clusters_of_h_m_genes=0;


if Q1_load_data
    load('YOUR_PATH\human_tese1and2_replicates_germcells.mat');
    load('YOUR_PATH\mouse_4m_replicates_germcells_ordered.mat');
    
    load('YOUR_PATH\rawdata\human_mouse_orthologue_pairs_20171018.mat');
    load('YOUR_PATH\h_m_stage_expression.mat');

    
        %setup colormap
    cd 'YOUR_PATH\cbrewer'
    %cbrewer()
    type_cm = cbrewer('qual','Paired',19);
    m_germtype_cm = cbrewer('qual','Paired',13);
    h_germtype_cm = cbrewer('qual','Paired',14);
    RdBu_cm = flipud(cbrewer('div','RdBu',50));
    RdBu_cm = RdBu_cm(5:45,:);
    genecluster_edgecm = [0.15 0.15 0.15;cbrewer('qual','Dark2',5)];
    genecluster_facecm = [0.5 0.5 0.5;cbrewer('qual','Set2',5)];
    Blues_cm = cbrewer('seq','Blues',50);
    BuPu_cm = cbrewer('seq','BuPu',50);
    PuBuGn_cm = cbrewer('seq','PuBuGn',50);
    YlGnBu_cm = cbrewer('seq','YlGnBu',50);
    cd 'YOUR_PATH'
    
        mouse_germcell_colormap = [    
    0.3281    0.1523    0.5313
    0.5000    0.4492    0.6719 
    0.6094    0.2188    0.3789
    0.7734    0.1758    0.2930
    0.9297    0.3906    0.2695 
    0.9844    0.6563    0.3672 
    0.9922    0.8750    0.5430 
    0.9922    0.9375    0.6523
    0.9336    0.9688    0.6484
    0.7305    0.8867    0.6250 
    0.4453    0.7734    0.6406
    0.3242    0.6719    0.6758
    0.1123    0.4668    0.7002];

    human_germcell_colormap = [    
    0.3281    0.1523    0.5313 %Sg-1
    0.5000    0.4492    0.6719 %Sg-2
    0.6094    0.2188    0.3789 %Sc-1
    0.7734    0.1758    0.2930 %Sc-2
    0.9297    0.3906    0.2695 %Sc-3
    0.9844    0.6563    0.3672 %Sc-4
    0.9922    0.8750    0.5430 %RS-1
    0.9922    0.9375    0.6523 %RS-2
    0.9336    0.9688    0.6484 %RS-3
    0.7305    0.8867    0.6250 %RS-4
    0.4453    0.7734    0.6406 %ES-1
    0.3242    0.6719    0.6758 %ES-2
    0.2109    0.5508    0.7266 %ES-3
    0.0137    0.3828    0.6738 %ES-4
];

end %Load data and colormaps



if Q2_conserve   
    %% Human_mouse conserved genes
    %define input files
    single_cells   =  1;
    single_stages  =  1;
    
    if single_cells
        H = h_germcell_sm_ordered;
        Hname = human_genename;
        Hstage = h_germ_IDX_ordered;
        M = m_germcell_sm_ordered;
        Mname = mouse_genename;
        Mstage = m_germ_IDX_ordered;
    end
    if single_stages
        H = h_cluster_expression;
        Hname = human_genename;
       % Hstage = [1: length(unique(h_germ_IDX_ordered))]';
        Hstage = [1: length(unique(h_germ_IDX))]';
        M = m_cluster_expression;
        Mname = mouse_genename;
       Mstage = [1: length(unique(m_germ_IDX_ordered))]';
       % Mstage = [1: length(unique(m_germ_IDX))]';
        
    end
    
    %setup the conserved gene list
    [~,i_H,~]=intersect(human_mouse_orthologue_pairs_20171018(:,1),Hname,'stable');
    [~,i_M,~]=intersect(human_mouse_orthologue_pairs_20171018(:,2),Mname,'stable');
    H_M_cons=intersect(i_H,i_M);
    [~,~,i_H]=intersect(human_mouse_orthologue_pairs_20171018(H_M_cons,1),Hname,'stable');
    [~,~,i_M]=intersect(human_mouse_orthologue_pairs_20171018(H_M_cons,2),Mname,'stable');

    
    %% plots of human-mouse conservation
    H_M_k3 = [H(i_H,:) M(i_M,:)];
    H_M_corr=corrcoef(F_normalize(F_normalize(H_M_k3,'median'),'FT'));
    figure;
    %mouse-human corrcoef heatmap of all single cells
    left = 0.15; bottom = 0.15; width = 0.75; height = 0.75;
    s1=subplot('position',[left bottom width height]); 
    imagesc(flipud(H_M_corr( (length(H(1,:))+1) :  (length(H(1,:)) + length(M(1,:)) ),1:length(H(1,:))  )),[0.3 0.8]);
    colormap(s1,RdBu_cm);
    grid off;
%    set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', 'xcolor', 'w', 'ycolor', 'w');
    set(gca,'xtick',''); set(gca,'ytick',''); 
    title('Human-mouse Correlation');hold on;
    %mouse stage bar
    left = 0.1; bottom = 0.15; width = 0.04; height = 0.75; 
    s=subplot('position',[left bottom width height]);
    imagesc(flipud(Mstage)); colormap(s,mouse_germcell_colormap);
    set(gca,'xtick',''); set(gca,'ytick','');hold on; ylabel('Mouse spermatogenesis');
    %human stage bar
    left = 0.15; bottom = 0.1; width = 0.75; height = 0.04;
    s=subplot('position',[left bottom width height]);
    imagesc(Hstage'); colormap(s,human_germcell_colormap);
    set(gca,'xtick',''); set(gca,'ytick','');xlabel('Human spermatogenesis');hold off;

    clear left right bottom width
    
    
    
end  %---relates to Fig.S2B



if Q3_joint_PCA
    %% define input files
    single_cells   =  1;
    single_stages  =  0;
    
    if single_cells
       % H = h_germcell_sm_ordered;
        H = h_germcell_sm;
        Hname = human_genename;
        Hstage = h_germ_IDX;
        M = m_germcell_sm_ordered;
        Mname = mouse_genename;
        Mstage = m_germ_IDX_ordered;
    end
    if single_stages
        H = h_cluster_expression;
        Hname = human_genename;
        Hstage = [1: length(unique(h_germ_IDX_ordered))]';
        M = m_cluster_expression;
        Mname = mouse_genename;
        Mstage = [1: length(unique(m_germ_IDX_ordered))]';
    end

    %setup the conserved gene list
    [~,i_H,~]=intersect(human_mouse_orthologue_pairs_20171018(:,1),Hname,'stable');
    [~,i_M,~]=intersect(human_mouse_orthologue_pairs_20171018(:,2),Mname,'stable');
    H_M_cons=intersect(i_H,i_M);
    [~,~,i_H]=intersect(human_mouse_orthologue_pairs_20171018(H_M_cons,1),Hname,'stable');
    [~,~,i_M]=intersect(human_mouse_orthologue_pairs_20171018(H_M_cons,2),Mname,'stable');
    
    %human-mouse conservation of all single cells
    H_M_k3 = [H(i_H,:) M(i_M,:)];
    H_M_conserbegene = human_mouse_orthologue_pairs_20171018(H_M_cons,:);
       h_idx = 1:length(H(1,:));
       m_idx = (length(H(1,:))+1) :  (length(H(1,:)) + length(M(1,:)));
   
    %% Joint PCA of all the human-mouse single cells
    
    A=F_normalize(H_M_k3,'tpm');
    A_log=F_normalize(H_M_k3,'tpm_log10');
    fano_factor  = var(A,0,2)./mean(A,2);
    PCA_input_gene = H_M_conserbegene(find(fano_factor>2*nanmean(fano_factor)),:);
    [hm_coeff, hm_score, hm_latent, hm_tsquared, hm_explained] = pca(A_log(find(fano_factor>2*nanmean(fano_factor)),:)');
    
        [~,xi]=sort(hm_coeff(:,2),'ascend');    
        PCA_input_gene(xi(1:20),2)
    
    %select the genes contributing to each species
    figure;
    left = 0.1; bottom = 0.55; width = 0.8; height = 0.4;
    f=subplot('position',[left bottom width height]); 
    [~,xi]=sort(hm_coeff(:,2),'ascend');    
    [~,~,xj]=intersect(PCA_input_gene(xi(1:20),1),H_M_conserbegene,'stable')
    imagesc(A_log(xj,:));colormap(f,RdBu_cm);
    yticks(1:length(xj));yticklabels(H_M_conserbegene(xj,1));xticks([]);
%     ax=gca
%     ax.YAxisLocation = 'right';
%     yticklabels(H_M_conserbegene(xj,2))
    left = 0.1; bottom = 0.12; width = 0.8; height = 0.4;
    f=subplot('position',[left bottom width height]); 
    [~,xi]=sort(hm_coeff(:,2),'descend');    
    [~,~,xj]=intersect(PCA_input_gene(xi(1:20),1),H_M_conserbegene,'stable')
    imagesc(A_log(xj,:));colormap(f,RdBu_cm);
    yticks(1:length(xj));yticklabels(H_M_conserbegene(xj));xticks([]);
    left = 0.1; bottom = 0.05; width = 0.8*(2094/3622); height = 0.05;
    f=subplot('position',[left bottom width height]); 
    imagesc(h_germ_IDX(h_germcell_pseudotimeorder,:)');
    colormap(f,human_germcell_colormap);yticks([]);xticks([]);
    left = 0.1+0.8*(2094/3622); bottom = 0.05; width = 0.8*(1528/3622); height = 0.05;
    f=subplot('position',[left bottom width height]); 
    imagesc(m_germ_IDX_ordered');colormap(f,mouse_germcell_colormap);yticks([]);xticks([]);
    
    colormap(human_germcell_colormap)
    
    %figure; histogram(hm_coeff(:,2),100)
    figure;    
    left = 0.1; bottom = 0.1; width = 0.55; height = 0.9;
    f=subplot('position',[left bottom width height]); 
    s=scatter3(hm_score(h_idx,1),hm_score(h_idx,2),hm_score(h_idx,3),30,h_pseudotime,'filled');
    s.MarkerEdgeColor='k';
    s.MarkerEdgeAlpha=0.2;
    s.MarkerFaceAlpha = 0.7; 
    hold on; 
    s=scatter3(hm_score(m_idx,1),hm_score(m_idx,2),hm_score(m_idx,3),30,m_pseudotime,'s','filled');
    xlabel({'PC1',hm_explained(1)}); ylabel({'PC2',hm_explained(2)});zlabel({'PC3',hm_explained(3)});
    colormap(YlGnBu_cm);
    s.MarkerEdgeColor='k';
    s.MarkerEdgeAlpha=0.2;
    s.MarkerFaceAlpha = 0.7;
    legend({'human','mouse'});hold off;
      %PC1-PC3 reflects speramtogenesis
        left = 0.7; bottom = 0.55; width = 0.25; height = 0.3;
        f=subplot('position',[left bottom width height]); 
        s=scatter(hm_score(h_idx,1),hm_score(h_idx,3),10,h_pseudotime,'filled');
        xlabel('PC1'); ylabel('PC3');
        colormap(YlGnBu_cm); hold on;  
        s.MarkerEdgeColor='k';
        s.MarkerEdgeAlpha=0.3;
        s.MarkerFaceAlpha = 1;      
        s=scatter(hm_score(m_idx,1),hm_score(m_idx,3),10,m_pseudotime,'s','filled');
                set(gca,'xtick',''); set(gca,'ytick','');
        s.MarkerEdgeColor='k';
        s.MarkerEdgeAlpha=0.3;
        s.MarkerFaceAlpha = 1; hold off;
          %PC1-PC2 reflects speramtogenesis
        left = 0.7; bottom = 0.15; width = 0.25; height = 0.3;
        f=subplot('position',[left bottom width height]); 
        s=scatter(hm_score(h_idx,1),hm_score(h_idx,2),10,h_pseudotime,'filled');
        xlabel('PC1'); ylabel('PC2');
        colormap(YlGnBu_cm); hold on; 
        s.MarkerEdgeColor='k';
        s.MarkerEdgeAlpha=0.3;
        s.MarkerFaceAlpha = 1;      
        s=scatter(hm_score(m_idx,1),hm_score(m_idx,2),10,m_pseudotime,'s','filled');
        set(gca,'xtick',''); set(gca,'ytick','');
        s.MarkerEdgeColor='k';
        s.MarkerEdgeAlpha=0.3;
        s.MarkerFaceAlpha = 1;  hold off;
    figure;
    subplot(1,2,1);
        s=scatter(hm_score(h_idx,1),hm_score(h_idx,3),10,h_pseudotime,'filled');
        xlabel('PC1'); ylabel('PC3');
        colormap(YlGnBu_cm); hold on; 
        s.MarkerEdgeColor='k';
        s.MarkerEdgeAlpha=0.3;
        s.MarkerFaceAlpha = 1;      
        s=scatter(hm_score(m_idx,1),hm_score(m_idx,3),10,m_pseudotime,'s','filled');
        s.MarkerEdgeColor='k';
        s.MarkerEdgeAlpha=0.3;
        s.MarkerFaceAlpha = 1;  hold off;
        set(gca,'color','none');
    subplot(1,2,2);
        s=scatter(hm_score(h_idx,1),hm_score(h_idx,2),10,h_pseudotime,'filled');
        xlabel('PC1'); ylabel('PC2');
        colormap(YlGnBu_cm); hold on; 
        s.MarkerEdgeColor='k';
        s.MarkerEdgeAlpha=0.3;
        s.MarkerFaceAlpha = 1;      
        s=scatter(hm_score(m_idx,1),hm_score(m_idx,2),10,m_pseudotime,'s','filled');
        s.MarkerEdgeColor='k';
        s.MarkerEdgeAlpha=0.3;
        s.MarkerFaceAlpha = 1;  hold off;
        set(gca,'color','none');legend({'human','mouse'});

end  %---relates to Fig.S2E-I

  

    cd 'YOUR_PATH\cbrewer'
    genecluster_edgecm = [0.15 0.15 0.15;cbrewer('qual','Dark2',5)];
    genecluster_facecm = [0.5 0.5 0.5;cbrewer('qual','Set2',5)];
    RdBu_cm = flipud(cbrewer('div','RdBu',50));
    RdBu_cm = RdBu_cm(5:45,:);
    cd 'YOUR_PATH\Analysis'



