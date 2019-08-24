function [matrix,column_names,row_names] = F_tsv_import(filename)
%%import filename.tsv as three matlab variables:
% filename should be like: 'filename.tsv'
% 1. matrix(gene*cell)
% 2. row_names
% 3. column_names

%%-->import all data
imported_data=importdata(filename);
matrix       = imported_data.data;
row_names    = imported_data.textdata(2:end,1);
column_names = imported_data.textdata(1,2:end);

