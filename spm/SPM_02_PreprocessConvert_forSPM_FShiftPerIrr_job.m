%-----------------------------------------------------------------------
% Job saved on 30-Mar-2026 15:47:04 by cfg_util (rev $Rev: 8183 $)
% spm SPM - SPM25 (25.01.02)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.meeg.convert.dataset = '<UNDEFINED>';
matlabbatch{1}.spm.meeg.convert.mode.continuous.readall = 1;
matlabbatch{1}.spm.meeg.convert.channels{1}.all = 'all';
matlabbatch{1}.spm.meeg.convert.outfile = '';
matlabbatch{1}.spm.meeg.convert.eventpadding = 0;
matlabbatch{1}.spm.meeg.convert.blocksize = 3276800;
matlabbatch{1}.spm.meeg.convert.checkboundary = 1;
matlabbatch{1}.spm.meeg.convert.saveorigheader = 0;
matlabbatch{1}.spm.meeg.convert.inputformat = 'autodetect';
matlabbatch{2}.spm.meeg.preproc.prepare.D(1) = cfg_dep('Conversion: Converted Datafile', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','Dfname'));
matlabbatch{2}.spm.meeg.preproc.prepare.task{1}.defaulteegsens.multimodal.nasfid = 'nas';
matlabbatch{2}.spm.meeg.preproc.prepare.task{1}.defaulteegsens.multimodal.lpafid = 'lpa';
matlabbatch{2}.spm.meeg.preproc.prepare.task{1}.defaulteegsens.multimodal.rpafid = 'rpa';
