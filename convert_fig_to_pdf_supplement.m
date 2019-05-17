% CONVERT_FIG_TO_PDF_SUPPLEMENT Convert FIG files to PDF files.
%
%   CONVERT_FIG_TO_PDF_SUPPLEMENT is a script we used to convert all the
%   figures in the supplementary material from FIG files to PDF files. It
%   requires the export_fig function [Al18], which is available on
%   MathWorks File Exchange.
%
% REFERENCES:
%
%   [Al18]  Y. Altman. export_fig, version 2.0.0.0, Available online, May
%           2018.
%           URL: https://www.mathworks.com/matlabcentral/fileexchange/23629-export_fig

fig_dir = 'figures/';
save_pdf_files = dir([fig_dir, 'S-experiment4*']);

for f = 1:length(save_pdf_files)
    uiopen([fig_dir, save_pdf_files(f).name], 1)
    pdf_file_name = [fig_dir, strrep(save_pdf_files(f).name, '.mat.fig', '')];
    export_fig(pdf_file_name, '-pdf', '-transparent');
    close all
end