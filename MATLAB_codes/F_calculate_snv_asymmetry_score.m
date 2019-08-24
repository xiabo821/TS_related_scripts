function [asymmetry_score,asym_score_type,snv_level,snv_type,output_genenames] = F_calculate_snv_asymmetry_score(...
                                             fn_snv_type_number,fn_base_frequency,...
                                             varargin)
%% [asymmetry_score,asym_score_type,snv_level,snv_type,output_genenames] = F_calculate_snv_asymmetry_score(
%                                                 fn_snv_number,...
%                                                 fn_base_frequency,...
%                                                 varargin);
%  F_compare_snv_asymmetry is wrote for 'transcriptional-scanning' manuscript.
%  Overall, this function use [SNV frequency] as input, [gene_cluster] as
%  grouping/staging method, [reference_group] as the reference for
%  comparison. Then this function plots the results of SNV frequency
%  between coding strand and template strand, as well as asymmetry scores,
%  and generates ranksum test p_values with Bonferroni correction.
%
%
%  Input:
%       fn_snv_type_number: file name of .tsv file (genes_by_SNVtypenumbers).   
%                           First column should be gene names. Other rows
%                           contains SNV number of each mutation type, including: 
%                           A-T,A-G,A-C,T-A,T-G,T-C,G-A,G-T,G-C,C-A,C-T,C-G
%       fn_base_frequency : file name of .tsv file which contains A/T/G/C 
%                           numbers of a given gene/sequence. First column 
%                           should be genenames. Note: base_frequency gene name  
%                           and snv_number gene names should be the same type. 
%  varargin:
%       'default_snv_asymmetry': SNV asymmetry calculation.
%                           default:{'AT','TA','AG','TC','TG','AC','CT','GA','GT','CA','CG','GC'}
%                           stands for: A>T/T>A; A>G/T>C; T>G/A>C; C>T/G>A; G>T/C>A; C>G/G>C;
%       'Subtraction'           : Option for involving substraction for the
%                                 fn_snv_type_number input.
%                                 Default is 'no'. 
%                                 Note: if option 'yes', do need to assign
%                                 fn_snv_type_number_ref and
%                                 fn_base_frequency_ref, which should have
%                                 exactly the same data organization with the
%                                 fn_snv_type_number file and
%                                 fn_base_frequency  file, respectively.
%       'fn_snv_type_number_ref': The optional input file used as reference
%                                 to subtract from the the fn_snv_type_number
%                                 input. Note this file should to have
%                                 exactly the same data organization with the
%                                 fn_snv_type_number file.
%       'fn_base_frequency_ref' : The optional input file used as reference
%                                 to subtract from the fn_base_frequency
%                                 input. Note this file should to have
%                                 exactly the same data organization with the
%                                 fn_base_frequency file.
%  
%  Output:
%       asymmetry_score     : gene-by-6type of SNV asymmetry score.
%                             Covering all genes.
%       asym_score_type     : corresponding to column names of snv_asymmetry_score.
%       snv_level           : gene-by-12type of SNV level, as defined by
%                             mutation numbers per gene per kb per 
%                             Covering all genes.
%       snv_type            : corresponding to column names of snv_level.
%       output_genenames    : corresponding to row names of asymmetry_score and snv_level.
%
%  Bo Xia
%


%% Parse defaul settings and setup inputs
fprintf('****************************************\n')
time1 = datetime;
fprintf('%s     Analysis started.\n',datestr(time1))

%default asymmetry calculations: 
%                       A>T/T>A;A>G/T>C;T>G/A>C;C>T/G>A;G>T/C>A;C>G/G>C;
default_snv_asymmetry = {'AT','TA','AG','TC','TG','AC','CT','GA','GT','CA','CG','GC'};

defaultSubtraction = 'no';
validSubtraction = {'no','yes'};
checkSubtraction = @(x) any(validatestring(x,validSubtraction));

defaultfn_snv_type_number_ref = {''};
defaultfn_base_frequency_ref = {''};

p = inputParser;
addOptional(p,'snv_asymmetry_order',default_snv_asymmetry)
addOptional(p,'Subtraction',defaultSubtraction,checkSubtraction)
addOptional(p,'fn_snv_type_number_ref',defaultfn_snv_type_number_ref)
addOptional(p,'fn_base_frequency_ref',defaultfn_base_frequency_ref)
parse(p,varargin{:})

%default numbers
default_snv_asymmetry = p.Results.snv_asymmetry_order;
defaultSubtraction = p.Results.Subtraction;
defaultfn_snv_type_number_ref = p.Results.fn_snv_type_number_ref;
defaultfn_base_frequency_ref = p.Results.fn_base_frequency_ref;

%setup input files
    %[SNV_number,SNV_genenames,SNV_type]
    [vnum,vgnames,vtype_names] = tsv_import(fn_snv_type_number);
    %[base_frequency,genenames,base_type]
    [vglength,vglenames,base_names] = tsv_import(fn_base_frequency);
    
    %input sample check if involves subtraction of input SNV_type numbers.
    if strcmp('yes',cellstr(defaultSubtraction))
        if strcmp('',cellstr(defaultfn_snv_type_number_ref)) | strcmp('',cellstr(defaultfn_base_frequency_ref)) 
            fprintf('Need to assign the "Subtraction" reference inputs.\n')
            return
        else
            [vnum_ref,vgnames_ref,vtype_names_ref] = tsv_import(defaultfn_snv_type_number_ref);
            [vglength_ref,vglenames_ref,base_names_ref] = tsv_import(defaultfn_base_frequency_ref);
            [~,a,b]=intersect(vgnames,vgnames_ref,'stable');
            [~,c,d]=intersect(vtype_names,vtype_names_ref,'stable');
            [~,e,f]=intersect(vglenames,vglenames_ref,'stable');
            [~,g,h]=intersect(base_names,base_names_ref,'stable');
            if sum(a-b)==0 & sum(c-d)==0 
                vnum     = vnum - vnum_ref;
            else
                fprintf('Inputs mismatch error! \n fn_snv_type_number and fn_snv_type_number_ref should have identical data organization.\n')
                return
            end
            if sum(e-f)==0 & sum(g-h)==0 
                vglength = vglength - vglength_ref;
            else
                fprintf('Inputs mismatch error! \n fn_base_frequency and fn_base_frequency_ref should have identical data organization.\n')
                return
            end
        end
    end

    %other default values
    asym_order  = default_snv_asymmetry;

%% set up SNV frequency (per kb per gene per SNV type)

[valid_genenames,i,j]=intersect(vgnames,vglenames,'stable');
vnum     = vnum(i,:);
vglength = vglength(j,:);
if length(j)==0
    fprintf('Error! Mismatched gene name of input SNV/base_frequency files!')
    return
end
%calculate SNV frequency per kilo-bases for each gene.
vnum_perkb = [];
vnum_perkb(:,strcmp('AT',asym_order)) = vnum(:,strcmp('AT',vtype_names))...
                                        ./   vglength(:,strcmp('A',base_names)) *1000;
vnum_perkb(:,strcmp('AG',asym_order)) = vnum(:,strcmp('AG',vtype_names))...
                                        ./   vglength(:,strcmp('A',base_names)) *1000;
vnum_perkb(:,strcmp('AC',asym_order)) = vnum(:,strcmp('AC',vtype_names))...
                                        ./   vglength(:,strcmp('A',base_names)) *1000;
vnum_perkb(:,strcmp('TA',asym_order)) = vnum(:,strcmp('TA',vtype_names))...
                                        ./   vglength(:,strcmp('T',base_names)) *1000;
vnum_perkb(:,strcmp('TG',asym_order)) = vnum(:,strcmp('TG',vtype_names))...
                                        ./   vglength(:,strcmp('T',base_names)) *1000;
vnum_perkb(:,strcmp('TC',asym_order)) = vnum(:,strcmp('TC',vtype_names))...
                                        ./   vglength(:,strcmp('T',base_names)) *1000;
vnum_perkb(:,strcmp('GA',asym_order)) = vnum(:,strcmp('GA',vtype_names))...
                                        ./   vglength(:,strcmp('G',base_names)) *1000;
vnum_perkb(:,strcmp('GT',asym_order)) = vnum(:,strcmp('GT',vtype_names))...
                                        ./   vglength(:,strcmp('G',base_names)) *1000;
vnum_perkb(:,strcmp('GC',asym_order)) = vnum(:,strcmp('GC',vtype_names))...
                                        ./   vglength(:,strcmp('G',base_names)) *1000;
vnum_perkb(:,strcmp('CA',asym_order)) = vnum(:,strcmp('CA',vtype_names))...
                                        ./   vglength(:,strcmp('C',base_names)) *1000;
vnum_perkb(:,strcmp('CT',asym_order)) = vnum(:,strcmp('CT',vtype_names))...
                                        ./   vglength(:,strcmp('C',base_names)) *1000;
vnum_perkb(:,strcmp('CG',asym_order)) = vnum(:,strcmp('CG',vtype_names))...
                                        ./   vglength(:,strcmp('C',base_names)) *1000;

valid_genenames = valid_genenames( ~any( isnan( vnum_perkb ) | isinf( vnum_perkb ), 2 ) );
vnum_perkb = vnum_perkb( ~any( isnan( vnum_perkb ) | isinf( vnum_perkb ), 2 ) ,: );


%% asymmetry score calculations
for i = 1:(length(asym_order)./2)
    asym_type(i) = strcat(asym_order(2*i-1),{' / '},asym_order(2*i));
    asymmetry_score(:,i) = log2((vnum_perkb(:,2*i-1)+1) ./ (vnum_perkb(:,2*i)+1) );
end

  
fprintf('     Finished.\n')
fprintf('****************************************\n')


%% outputs
asymmetry_score     = asymmetry_score;
asym_score_type     = asym_type;
snv_level           = vnum_perkb;
snv_type            = asym_order;
output_genenames    = valid_genenames;




%% functions
function [matrix,row_names,column_names] = tsv_import(filename)
%%import filename.tsv as three matlab variables:
% filename should be like: 'filename.tsv'
% 1. matrix(gene*cell)
% 2. row_names
% 3. column_names

%%-->import all data
imported_data= importdata(filename);
matrix       = imported_data.data;
row_names    = imported_data.textdata(2:end,1);
column_names = imported_data.textdata(1,2:end);