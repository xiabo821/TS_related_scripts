function [filtered_matrix,filtered_cellnames,gene_names] = F_inDrop_filtercells(inDrop_raw_input_filename,varargin)
%% [filtered_matrix,filtered_cellnames,gene_names] = F_inDrop_filter(inDrop_raw_input_filename,...
%                                                                    varargin)
%   e.g.:
%   [clean_max,cellnames,genenames] = F_inDrop_filter(inDrop_raw_input_filename,...
%                                                     'UMI_cutoff',0,...
%                                                     'Gene_number_cutoff',1000,...
%                                                     'Mito_ratio_cutoff',0.2,...
%                                                     'Ribo_ratio_cutoff',0.2,...
%                                                     'Remove_doublet','yes',...
%                                                     );
%   
%   This function work with the inDrop raw input files, or any expression
%   matrix recognizable by importdata.m function.
%   The raw input file usually is tab delimited .tsv/.txt gene expression file
%   from inDrop pipeline or any expression matrix with fist column being
%   gene names and first row being the cell names.
%   This function firstly provides some general states of the raw results.
%   Then the function removes cells with bad quality and/or being potential doublets.
%
%  Required functions: 
%           importdata.m            Copyright@ MATLAB
%
%  INPUT:
%  inDrop_raw_input_filename: file path, better full path. 
%                             Usually is tab delimited .tsv/.txt gene expression file
%
%  varargin:
%     VARARGIN, DEFAUT >>> INFO
%     'UMI_cutoff',0 >>> detect more than 0 UMIs; (no UMI requirement)
%                  Option: suggest 1000 UMIs, or any integer.
%     'Gene_number_cutoff',1000 >>> detect more than 1000 genes;
%                  Option: suggest 1000 genes, or any integer.
%     'Mito_ratio_cutoff',0.2 >>> detect less than 20% mitochondrial genes 
%                  Option: suggest 0.2 (20%) genes, or any number in [0 1].
%     'Ribo_ratio_cutoff',0.2 >>> detect less than 20% ribosomal genes
%                  Option: suggest 0.2 (20%) genes, or any number in [0 1].
%     'Remove_doublet','yes' >>> filter the cells which is outlier 
%                                (MATLAB 'isoutlier' function with 'mean' method)
%                                in the highly-transcripted cells.
%                                Highly-transcripted cells are defined as the 
%                                4th quarter of total transcripts of each cell 
%                                (after gene number and mito/ribo filtering). 
%                  Option: 'no'
%   
%   Note: this function has only been tested on human and mouse samples.
%   Be careful if applying to other species, especially mitochondrial and
%   ribosomal genes.
%   
%   Bo Xia (since 08/21/2018)

%% import data
time1 = datetime;
fprintf('*********************************************************\n') 
fprintf('Your started at %s.\n',string(time1));

fn = inDrop_raw_input_filename;

    a          = importdata(fn);
    exp_matrix = a.data;
    gene_names = a.textdata(2:end,1);
    cell_names = a.textdata(1,2:end);
    
fprintf('Data imported.\n');
fprintf('\n')

%% set up varargins
%default values:
umi_cutoff         = 0; %should detect more than UMI_number_cutoff UMIs
gene_number_cutoff = 1000; %should detect more than gene_number_cutoff genes
mito_ratio_cutoff  = 0.2; % should have less that mito_ratio_cutoff mito_genes
ribo_ratio_cutoff  = 0.2; % should have less that ribo_ratio_cutoff ribo_genes
default_remove_doublet = 'yes';
validOption = {'no','yes'};
checkOption = @(x) any(validatestring(x,validOption));

p = inputParser;
addOptional(p,'UMI_cutoff',umi_cutoff,@isnumeric)
addOptional(p,'Gene_number_cutoff',gene_number_cutoff,@isnumeric)
addOptional(p,'Mito_ratio_cutoff',mito_ratio_cutoff,@isnumeric)
addOptional(p,'Ribo_ratio_cutoff',ribo_ratio_cutoff,@isnumeric)
addOptional(p,'Remove_doublet',default_remove_doublet,checkOption)
parse(p,varargin{:})

%default values
umi_cutoff         = p.Results.UMI_cutoff;
gene_number_cutoff = p.Results.Gene_number_cutoff;
mito_ratio_cutoff  = p.Results.Mito_ratio_cutoff;
ribo_ratio_cutoff  = p.Results.Ribo_ratio_cutoff;
default_remove_doublet  = p.Results.Remove_doublet;


%% plot basic numbers of transcripts
    A = exp_matrix;
    figure;
    subplot(2,1,1); histogram(log10(sum(A)),100,'EdgeAlpha',0.5,'FaceColor',[0 0.5 0.5]);
    ylabel('Cell numbers')
    xlabel('Total transcript numbers (log10)')
    A(A~=0) = 1;
    subplot(2,1,2); histogram(sum(A),100,'EdgeAlpha',0.5,'FaceColor',[0.75 0.68 0.84]);
    ylabel('Cell numbers')
    xlabel('Detected gene numbers') 
    suptitle('Transcripts numbers (top) and Gene numbers (bottom) of all cells');
fprintf('Check out the plots of raw transcripts information of all cells.\n');
fprintf('\n')

    
%% plot mitochondial/ribosomal genes ratio

    A = exp_matrix;
    %mito_ratio. Work for human and mouse
    mito = [strmatch('MT-',gene_names);strmatch('mt-',gene_names)];
    mito_ratio=(sum(A(mito,:))./sum(A));
    %ribo_ratio. Work for human and mouse
    ribo = [strmatch('RPL',gene_names);strmatch('RPS',gene_names);...
        strmatch('Rpl',gene_names);strmatch('Rps',gene_names)]; 
    ribo_ratio=(sum(A(ribo,:))./sum(A));
    
    %plot mitochondrial genes and ribosomal genes
    figure;
    %mitochondrial transcripts ratio
    subplot(2,2,1); histogram(mito_ratio,50,'Orientation', 'horizontal',...
        'EdgeAlpha',0.5,'FaceColor',[0.89 0.10 0.11]);
    xlabel('Cell numbers') 
    ylabel('Mito-transcripts ratio')
    %mitochondrial transcripts vs total transcripts
    subplot(2,2,2); scatter(log10(sum(A)),mito_ratio,5,'MarkerEdgeColor','none',...
        'MarkerFaceAlpha',0.5,'MarkerFaceColor',[0.89 0.10 0.11]);
    xlabel('Total transcripts per cell(log10)') 
    ylabel('Mito-transcripts ratio')
    %Ribosomal transcripts ratio 
    subplot(2,2,3); histogram(ribo_ratio,50,'Orientation', 'horizontal',...
        'EdgeAlpha',0.5,'FaceColor',[0.122 0.472 0.709]);
    xlabel('Cell numbers') 
    ylabel('Ribo-transcripts ratio')
    %Ribosomal transcripts ratio
    subplot(2,2,4); scatter(log10(sum(A)),ribo_ratio,5,'MarkerEdgeColor','none',...
        'MarkerFaceAlpha',0.5,'MarkerFaceColor',[0.122 0.472 0.709]);
    xlabel('Total transcripts per cell (log10)') 
    ylabel('Ribo-transcripts ratio')
    suptitle('Mitochondrial/ribosomal transcripts information of all cells')
fprintf('Check out the plots of mitochondrial/ribosomal transcripts information.\n');
fprintf('\n')

    
%% Filter cells 
fprintf('Start to filter cells.\n');

    %Display your QC/filtering criteria
    fprintf('\t Your QC/filtering criteria:\n');
    fprintf('\t\t UMI_cutoff:         %d\n',umi_cutoff);
    fprintf('\t\t Gene_number_cutoff: %d\n',gene_number_cutoff);
    fprintf('\t\t Mito_ratio_cutoff:  %d\n',mito_ratio_cutoff);
    fprintf('\t\t Ribo_ratio_cutoff:  %d\n',ribo_ratio_cutoff);
    fprintf('\t\t Remove_doublet:     %s\n',default_remove_doublet);

    % Filter by detected gene number and mitochondial/ribosomal genes ratios
    A = exp_matrix;
    A(A~=0) = 1;

    good_quality_cell = intersect( intersect(find(mito_ratio<mito_ratio_cutoff),...
                                     find(ribo_ratio<ribo_ratio_cutoff) ),...
                                   intersect(find(sum(exp_matrix)>UMI_cutoff),...
                                     find(sum(A)>gene_number_cutoff))     );
fprintf('Input cell number: %d .\n',length(mito_ratio)) 
fprintf('Your cutoff:\n')
fprintf('   detects more than %d transcripts(UMI);\n',UMI_cutoff)
fprintf('   express more than %d genes;\n',gene_number_cutoff)
fprintf('   mitochondial transcriptsless < %d percents;\n',...
                                           mito_ratio_cutoff*100)
fprintf('   ribosomal gene transcripts < %d percents.\n',...
                                           ribo_ratio_cutoff*100)
fprintf('You got %d cells with good sequencing quality.\n',length(good_quality_cell)) 
fprintf('\n')
fprintf('\n')
    
    % Filter potential doublelets
    fprintf('Out of the %d cells with good quality,\n',length(good_quality_cell)) 

    if strcmp('yes',cellstr(default_remove_doublet))
    A = exp_matrix(:,good_quality_cell);
    A_sum = sum(A);
%     figure; histogram(A_sum,100,...
%         'EdgeAlpha',0.5,'FaceColor',[0.122 0.472 0.709]);
    %filter the cells which is outlier in the highly-transcripted cells.
    %highly-transcripted cells are defined as the 4th quarter of total 
    %transcripts of each cell (after gene number and mito/ribo filtering).
    Fourth_quarter_cells = find(A_sum >= quantile(A_sum,0.75));
    possible_doublets  = intersect( find(isoutlier(A_sum(Fourth_quarter_cells),'mean')),...
                                    find(A_sum(Fourth_quarter_cells) > mean(A_sum(Fourth_quarter_cells)))...
                                    );
    output_cell = good_quality_cell ( ...
                     setdiff(1:length(good_quality_cell),possible_doublets) );
    fprintf('%d cells is removed as potential doublets.\n',length(possible_doublets)) 
    fprintf('Finally, you have %d cells.\n',(length(good_quality_cell)-length(possible_doublets))) 
    else
        output_cell = good_quality_cell;
        fprintf(' without performing potential doublets checking, \n',length(possible_doublets)) 
        fprintf('Finally, you have %d cells.\n',length(output_cell)) 
    end
    

fprintf('\n')

                    
%% Export output files
    filtered_matrix    = exp_matrix(:,output_cell);
    gene_names         = a.textdata(2:end,1);
    filtered_cellnames = cell_names(output_cell);

%% Plot the states of the output cells    
fprintf('Plot some of the stats of good cells.\n');
    A = filtered_matrix;

    figure;
    subplot(2,1,1); 
    mean_transcript = round(mean(sum(A)));
    median_transcript = round(median(sum(A)));    
    histogram(log10(sum(A)),50,...
        'EdgeAlpha',0.5,'FaceColor',[0 0.5 0.5]);hold on
    line([log10(mean_transcript) log10(mean_transcript)],...
         [0 max(get(gca,'YLim'))],...
         'Color',[227 26 28]/256,'LineStyle','--');
        str = ['\leftarrow mean=',num2str(mean_transcript)];
        text(log10(mean_transcript),max(get(gca,'YLim'))*0.8, str)
    line([log10(median_transcript) log10(median_transcript)],...
         [0 max(get(gca,'YLim'))],...
         'Color',[227 26 28]/256,'LineStyle','--');
        str = ['\leftarrow median=',num2str(median_transcript)];
        text(log10(median_transcript),max(get(gca,'YLim'))*0.9, str)
    hold off
    xlabel('Total transcripts (log10)') 
    ylabel('Cell numbers')

    A(A~=0) = 1;
    mean_transcript = round(mean(sum(A)));
    median_transcript = round(median(sum(A)));    
    subplot(2,1,2); histogram(sum(A),50,...
        'EdgeAlpha',0.5,'FaceColor',[0.75 0.68 0.84]);
    line([mean_transcript mean_transcript],...
         [0 max(get(gca,'YLim'))],...
         'Color',[227 26 28]/256,'LineStyle','--');
        str = ['\leftarrow mean=',num2str(mean_transcript)];
        text(mean_transcript,max(get(gca,'YLim'))*0.8, str)
    line([median_transcript median_transcript],...
         [0 max(get(gca,'YLim'))],...
         'Color',[227 26 28]/256,'LineStyle','--');
        str = ['\leftarrow median=',num2str(median_transcript)];
        text(median_transcript,max(get(gca,'YLim'))*0.9, str)
    hold off
    xlabel('Detected gene numbers') 
    ylabel('Cell numbers')
    txt = ['of the ',num2str( length(filtered_matrix(1,:)) ),' cells passing filter'];
    suptitle({'Transcripts numbers (log10, top) and Gene numbers (bottom)',...
              txt});
    
fprintf('You are done! In total it takes: %s.\n',string(datetime - time1));
fprintf('\n')
fprintf('\n')


