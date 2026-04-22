%-----------------------------------------------------------------------
% shortened contrast manager script
%-----------------------------------------------------------------------
F.SPM_File='\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spm_ANOVA_FlexibleFactorial_onplyprecue_attaveraged_only26Hz\SPM.mat';
F.SPM_File='\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spm_ANOVA_FlexibleFactorial_onplyprecue_attaveraged\SPM.mat';

matlabbatch = [];
matlabbatch{1}.spm.stats.con.spmmat = {F.SPM_File};

matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Identity';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.fcon.name = 'MAIN LUMINANCE';
matlabbatch{1}.spm.stats.con.consess{2}.fcon.weights = [1 -1];
matlabbatch{1}.spm.stats.con.consess{2}.fcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Offset > Isolum';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [1 -1];
matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'Isolum > Offset';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [-1 1];
matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.delete = 1;
spm_jobman('run',matlabbatch);