%-----------------------------------------------------------------------
% shortened contrast manager script
%-----------------------------------------------------------------------
F.SPM_File='\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spm_ANOVA_FlexibleFactorial_onplypostcue\SPM.mat';
F.SPM_File='\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spm_ANOVA_FlexibleFactorial_onplypostcue_allssveps\SPM.mat';
% F.SPM_File='\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spm_ANOVA_FlexibleFactorial_onplypostcue_only26Hz\SPM.mat';
% F.SPM_File='\\smbone.dom.uni-leipzig.de\FFL\AllgPsy\experimental_data\2024_FShiftPerIrr\eeg\SPM\spm_ANOVA_FlexibleFactorial_onplyprecue_only26Hz\SPM.mat';

matlabbatch = [];
matlabbatch{1}.spm.stats.con.spmmat = {F.SPM_File};

matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Identity';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 1 1 1];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.fcon.name = 'MAIN Att vs Unatt';
matlabbatch{1}.spm.stats.con.consess{2}.fcon.weights = [1 1 -1 -1];
matlabbatch{1}.spm.stats.con.consess{2}.fcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{3}.fcon.name = 'MAIN Isolum vs Offset';
matlabbatch{1}.spm.stats.con.consess{3}.fcon.weights = [-1 1 -1 1];
matlabbatch{1}.spm.stats.con.consess{3}.fcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{4}.fcon.name = 'INTERACTION (Isol_Att>Isol_Unatt) vs (Off_Att>Off_Unatt)';
matlabbatch{1}.spm.stats.con.consess{4}.fcon.weights = [-1 1 1 -1];
matlabbatch{1}.spm.stats.con.consess{4}.fcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'Offset > Isolum';
matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [1 -1 1 -1];
matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'Isolum > Offset';
matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [-1 1 -1 1];
matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'Offset: Att > Unatt';
matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [1 0 -1 0];
matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'Isolum: Att > Unatt';
matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = [0 1 0 -1];
matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{9}.tcon.name = 'Att: Isolum > Offset';
matlabbatch{1}.spm.stats.con.consess{9}.tcon.weights = [-1 1 0 0];
matlabbatch{1}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{10}.tcon.name = '(Isol_Att>Isol_Unatt) > (Off_Att>Off_Unatt)';
matlabbatch{1}.spm.stats.con.consess{10}.tcon.weights = [-1 1 1 -1];
matlabbatch{1}.spm.stats.con.consess{10}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.delete = 1;
spm_jobman('run',matlabbatch);