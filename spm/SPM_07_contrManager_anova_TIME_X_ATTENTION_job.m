%-----------------------------------------------------------------------
% shortened contrast manager script
%-----------------------------------------------------------------------
F.SPM_File='\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spm_ANOVA_FlexibleFactorial_OnlyIsolum\SPM.mat';

matlabbatch = [];
matlabbatch{1}.spm.stats.con.spmmat = {F.SPM_File};

matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Identity';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 1 1 1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.fcon.name = 'MAIN Att vs Unatt';
matlabbatch{1}.spm.stats.con.consess{2}.fcon.weights = [1 1 -1 -1];
matlabbatch{1}.spm.stats.con.consess{2}.fcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{3}.fcon.name = 'MAIN Early vs Late';
matlabbatch{1}.spm.stats.con.consess{3}.fcon.weights = [-1 1 -1 1];
matlabbatch{1}.spm.stats.con.consess{3}.fcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{4}.fcon.name = 'INTERACTION (Early_Att>Early_Unatt) vs (Late_Att>Late_Unatt)';
matlabbatch{1}.spm.stats.con.consess{4}.fcon.weights = [-1 1 1 -1];
matlabbatch{1}.spm.stats.con.consess{4}.fcon.sessrep = 'none';


matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'Early: Att > Unatt';
matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [1 0 -1 0];
matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'Late: Att > Unatt';
matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [0 1 0 -1];
matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'INTERACTION (Late_Att>Late_Unatt) > (Early_Att>Early_Unatt)';
matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [-1 1 1 -1];
matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.delete = 1;

spm_jobman('run',matlabbatch);