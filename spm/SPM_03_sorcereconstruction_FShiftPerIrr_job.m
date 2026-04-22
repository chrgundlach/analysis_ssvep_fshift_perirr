%-----------------------------------------------------------------------
% Job saved on 31-Mar-2026 13:42:16 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM25 (25.01.02)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.meeg.source.headmodel.D = {'<UNDEFINED>'};
matlabbatch{1}.spm.meeg.source.headmodel.val = 1;
matlabbatch{1}.spm.meeg.source.headmodel.comment = '';
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshes.template = 1;
matlabbatch{1}.spm.meeg.source.headmodel.meshing.meshres = 2;
matlabbatch{1}.spm.meeg.source.headmodel.coregistration.coregdefault = 1;
matlabbatch{1}.spm.meeg.source.headmodel.forward.eeg = 'EEG BEM';
matlabbatch{1}.spm.meeg.source.headmodel.forward.meg = 'Single Shell';
matlabbatch{2}.spm.meeg.source.invert.D(1) = cfg_dep('Head model specification: M/EEG dataset(s) with a forward model', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','D'));
matlabbatch{2}.spm.meeg.source.invert.val = 1;
matlabbatch{2}.spm.meeg.source.invert.whatconditions.all = 1;
matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.invtype = 'GS';
matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.woi = [-Inf Inf];
matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.foi = [0 256];
matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.hanning = 1;
matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.priors.priorsmask = {''};
matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.priors.space = 1;
matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.restrict.locs = zeros(0, 3);
matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.restrict.radius = 32;
matlabbatch{2}.spm.meeg.source.invert.isstandard.custom.restrict.mask = {''};
matlabbatch{2}.spm.meeg.source.invert.modality = {'All'};
