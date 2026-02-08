%
% Status : main Dynare file
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

tic0 = tic;
% Define global variables.
global M_ options_ oo_ estim_params_ bayestopt_ dataset_ dataset_info estimation_info
options_ = [];
M_.fname = 'HPJedit_baseline_RBC';
M_.dynare_version = '6.5';
oo_.dynare_version = '6.5';
options_.dynare_version = '6.5';
%
% Some global variables initialization
%
global_initialization;
M_.exo_names = cell(1,1);
M_.exo_names_tex = cell(1,1);
M_.exo_names_long = cell(1,1);
M_.exo_names(1) = {'T_shock'};
M_.exo_names_tex(1) = {'T\_shock'};
M_.exo_names_long(1) = {'T_shock'};
M_.endo_names = cell(10,1);
M_.endo_names_tex = cell(10,1);
M_.endo_names_long = cell(10,1);
M_.endo_names(1) = {'Z'};
M_.endo_names_tex(1) = {'Z'};
M_.endo_names_long(1) = {'Z'};
M_.endo_names(2) = {'R'};
M_.endo_names_tex(2) = {'R'};
M_.endo_names_long(2) = {'R'};
M_.endo_names(3) = {'W'};
M_.endo_names_tex(3) = {'W'};
M_.endo_names_long(3) = {'W'};
M_.endo_names(4) = {'L'};
M_.endo_names_tex(4) = {'L'};
M_.endo_names_long(4) = {'L'};
M_.endo_names(5) = {'K'};
M_.endo_names_tex(5) = {'K'};
M_.endo_names_long(5) = {'K'};
M_.endo_names(6) = {'I'};
M_.endo_names_tex(6) = {'I'};
M_.endo_names_long(6) = {'I'};
M_.endo_names(7) = {'Y'};
M_.endo_names_tex(7) = {'Y'};
M_.endo_names_long(7) = {'Y'};
M_.endo_names(8) = {'C'};
M_.endo_names_tex(8) = {'C'};
M_.endo_names_long(8) = {'C'};
M_.endo_names(9) = {'G'};
M_.endo_names_tex(9) = {'G'};
M_.endo_names_long(9) = {'G'};
M_.endo_names(10) = {'T'};
M_.endo_names_tex(10) = {'T'};
M_.endo_names_long(10) = {'T'};
M_.endo_partitions = struct();
M_.param_names = cell(16,1);
M_.param_names_tex = cell(16,1);
M_.param_names_long = cell(16,1);
M_.param_names(1) = {'Phi'};
M_.param_names_tex(1) = {'Phi'};
M_.param_names_long(1) = {'Phi'};
M_.param_names(2) = {'delta'};
M_.param_names_tex(2) = {'delta'};
M_.param_names_long(2) = {'delta'};
M_.param_names(3) = {'sigma'};
M_.param_names_tex(3) = {'sigma'};
M_.param_names_long(3) = {'sigma'};
M_.param_names(4) = {'alpha'};
M_.param_names_tex(4) = {'alpha'};
M_.param_names_long(4) = {'alpha'};
M_.param_names(5) = {'rho'};
M_.param_names_tex(5) = {'rho'};
M_.param_names_long(5) = {'rho'};
M_.param_names(6) = {'chi'};
M_.param_names_tex(6) = {'chi'};
M_.param_names_long(6) = {'chi'};
M_.param_names(7) = {'Z_ss'};
M_.param_names_tex(7) = {'Z\_ss'};
M_.param_names_long(7) = {'Z_ss'};
M_.param_names(8) = {'R_ss'};
M_.param_names_tex(8) = {'R\_ss'};
M_.param_names_long(8) = {'R_ss'};
M_.param_names(9) = {'W_ss'};
M_.param_names_tex(9) = {'W\_ss'};
M_.param_names_long(9) = {'W_ss'};
M_.param_names(10) = {'L_ss'};
M_.param_names_tex(10) = {'L\_ss'};
M_.param_names_long(10) = {'L_ss'};
M_.param_names(11) = {'K_ss'};
M_.param_names_tex(11) = {'K\_ss'};
M_.param_names_long(11) = {'K_ss'};
M_.param_names(12) = {'I_ss'};
M_.param_names_tex(12) = {'I\_ss'};
M_.param_names_long(12) = {'I_ss'};
M_.param_names(13) = {'Y_ss'};
M_.param_names_tex(13) = {'Y\_ss'};
M_.param_names_long(13) = {'Y_ss'};
M_.param_names(14) = {'C_ss'};
M_.param_names_tex(14) = {'C\_ss'};
M_.param_names_long(14) = {'C_ss'};
M_.param_names(15) = {'G_ss'};
M_.param_names_tex(15) = {'G\_ss'};
M_.param_names_long(15) = {'G_ss'};
M_.param_names(16) = {'T_ss'};
M_.param_names_tex(16) = {'T\_ss'};
M_.param_names_long(16) = {'T_ss'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 10;
M_.param_nbr = 16;
M_.orig_endo_nbr = 10;
M_.aux_vars = [];
M_.Sigma_e = zeros(1, 1);
M_.Correlation_matrix = eye(1, 1);
M_.H = 0;
M_.Correlation_matrix_ME = 1;
M_.sigma_e_is_diagonal = true;
M_.det_shocks = [];
M_.surprise_shocks = [];
M_.learnt_shocks = [];
M_.learnt_endval = [];
M_.heteroskedastic_shocks.Qvalue_orig = [];
M_.heteroskedastic_shocks.Qscale_orig = [];
M_.matched_irfs = {};
M_.matched_irfs_weights = {};
options_.linear = false;
options_.block = false;
options_.bytecode = false;
options_.use_dll = false;
options_.ramsey_policy = false;
options_.discretionary_policy = false;
M_.eq_nbr = 10;
M_.ramsey_orig_eq_nbr = 0;
M_.ramsey_orig_endo_nbr = 0;
M_.set_auxiliary_variables = exist(['./+' M_.fname '/set_auxiliary_variables.m'], 'file') == 2;
M_.epilogue_names = {};
M_.epilogue_var_list_ = {};
M_.orig_maximum_endo_lag = 1;
M_.orig_maximum_endo_lead = 1;
M_.orig_maximum_exo_lag = 0;
M_.orig_maximum_exo_lead = 0;
M_.orig_maximum_exo_det_lag = 0;
M_.orig_maximum_exo_det_lead = 0;
M_.orig_maximum_lag = 1;
M_.orig_maximum_lead = 1;
M_.orig_maximum_lag_with_diffs_expanded = 1;
M_.lead_lag_incidence = [
 1 3 0;
 0 4 13;
 0 5 0;
 0 6 0;
 2 7 0;
 0 8 0;
 0 9 0;
 0 10 14;
 0 11 0;
 0 12 15;]';
M_.nstatic = 5;
M_.nfwrd   = 3;
M_.npred   = 2;
M_.nboth   = 0;
M_.nsfwrd   = 3;
M_.nspred   = 2;
M_.ndynamic   = 5;
M_.dynamic_tmp_nbr = [3; 0; 0; 0; ];
M_.equations_tags = {
  1 , 'name' , '1' ;
  2 , 'name' , '2' ;
  3 , 'name' , 'Y' ;
  4 , 'name' , 'W' ;
  5 , 'name' , 'R' ;
  6 , 'name' , 'K' ;
  7 , 'name' , '7' ;
  8 , 'name' , 'G' ;
  9 , 'name' , 'Z' ;
  10 , 'name' , 'T' ;
};
M_.mapping.Z.eqidx = [3 9 ];
M_.mapping.R.eqidx = [1 5 ];
M_.mapping.W.eqidx = [2 4 ];
M_.mapping.L.eqidx = [2 3 4 ];
M_.mapping.K.eqidx = [3 5 6 ];
M_.mapping.I.eqidx = [6 7 ];
M_.mapping.Y.eqidx = [3 4 5 7 ];
M_.mapping.C.eqidx = [1 2 7 8 ];
M_.mapping.G.eqidx = [7 8 ];
M_.mapping.T.eqidx = [1 2 8 10 ];
M_.mapping.T_shock.eqidx = [10 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.block_structure.time_recursive = false;
M_.block_structure.block(1).Simulation_Type = 1;
M_.block_structure.block(1).endo_nbr = 2;
M_.block_structure.block(1).mfs = 2;
M_.block_structure.block(1).equation = [ 9 10];
M_.block_structure.block(1).variable = [ 1 10];
M_.block_structure.block(1).is_linear = true;
M_.block_structure.block(1).NNZDerivatives = 3;
M_.block_structure.block(1).bytecode_jacob_cols_to_sparse = [1 3 4 ];
M_.block_structure.block(2).Simulation_Type = 8;
M_.block_structure.block(2).endo_nbr = 8;
M_.block_structure.block(2).mfs = 6;
M_.block_structure.block(2).equation = [ 4 8 3 7 2 6 5 1];
M_.block_structure.block(2).variable = [ 3 9 7 6 4 5 2 8];
M_.block_structure.block(2).is_linear = false;
M_.block_structure.block(2).NNZDerivatives = 20;
M_.block_structure.block(2).bytecode_jacob_cols_to_sparse = [4 0 0 7 8 9 10 11 12 17 18 ];
M_.block_structure.block(1).g1_sparse_rowval = int32([]);
M_.block_structure.block(1).g1_sparse_colval = int32([]);
M_.block_structure.block(1).g1_sparse_colptr = int32([]);
M_.block_structure.block(2).g1_sparse_rowval = int32([1 4 5 1 2 3 5 2 4 1 3 4 5 2 3 6 6 6 ]);
M_.block_structure.block(2).g1_sparse_colval = int32([4 4 4 7 7 7 7 8 8 9 9 10 11 12 12 12 17 18 ]);
M_.block_structure.block(2).g1_sparse_colptr = int32([1 1 1 1 4 4 4 8 10 12 13 14 17 17 17 17 17 18 19 ]);
M_.block_structure.variable_reordered = [ 1 10 3 9 7 6 4 5 2 8];
M_.block_structure.equation_reordered = [ 9 10 4 8 3 7 2 6 5 1];
M_.block_structure.incidence(1).lead_lag = -1;
M_.block_structure.incidence(1).sparse_IM = [
 3 5;
 5 5;
 6 5;
 9 1;
];
M_.block_structure.incidence(2).lead_lag = 0;
M_.block_structure.incidence(2).sparse_IM = [
 1 8;
 1 10;
 2 3;
 2 4;
 2 8;
 2 10;
 3 1;
 3 4;
 3 7;
 4 3;
 4 4;
 4 7;
 5 2;
 5 7;
 6 5;
 6 6;
 7 6;
 7 7;
 7 8;
 7 9;
 8 8;
 8 9;
 8 10;
 9 1;
 10 10;
];
M_.block_structure.incidence(3).lead_lag = 1;
M_.block_structure.incidence(3).sparse_IM = [
 1 2;
 1 8;
 1 10;
];
M_.block_structure.dyn_tmp_nbr = 2;
M_.state_var = [1 5 ];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(10, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(16, 1);
M_.endo_trends = struct('deflator', cell(10, 1), 'log_deflator', cell(10, 1), 'growth_factor', cell(10, 1), 'log_growth_factor', cell(10, 1));
M_.NNZDerivatives = [33; -1; -1; ];
M_.dynamic_g1_sparse_rowval = int32([9 3 5 6 3 9 5 2 4 2 3 4 6 6 7 3 4 5 7 1 2 7 8 7 8 1 2 8 10 1 1 1 10 ]);
M_.dynamic_g1_sparse_colval = int32([1 5 5 5 11 11 12 13 13 14 14 14 15 16 16 17 17 17 17 18 18 18 18 19 19 20 20 20 20 22 28 30 31 ]);
M_.dynamic_g1_sparse_colptr = int32([1 2 2 2 2 5 5 5 5 5 5 7 8 10 13 14 16 20 24 26 30 30 31 31 31 31 31 31 32 32 33 34 ]);
M_.lhs = {
'1/((1+T)*C)'; 
'chi/(1-L)'; 
'Y'; 
'W'; 
'R'; 
'K'; 
'Y'; 
'G'; 
'Z'; 
'T'; 
};
M_.static_tmp_nbr = [3; 0; 0; 0; ];
M_.block_structure_stat.block(1).Simulation_Type = 3;
M_.block_structure_stat.block(1).endo_nbr = 1;
M_.block_structure_stat.block(1).mfs = 1;
M_.block_structure_stat.block(1).equation = [ 9];
M_.block_structure_stat.block(1).variable = [ 1];
M_.block_structure_stat.block(2).Simulation_Type = 1;
M_.block_structure_stat.block(2).endo_nbr = 1;
M_.block_structure_stat.block(2).mfs = 1;
M_.block_structure_stat.block(2).equation = [ 10];
M_.block_structure_stat.block(2).variable = [ 10];
M_.block_structure_stat.block(3).Simulation_Type = 6;
M_.block_structure_stat.block(3).endo_nbr = 8;
M_.block_structure_stat.block(3).mfs = 8;
M_.block_structure_stat.block(3).equation = [ 3 4 5 6 7 8 1 2];
M_.block_structure_stat.block(3).variable = [ 4 3 5 6 7 9 2 8];
M_.block_structure_stat.variable_reordered = [ 1 10 4 3 5 6 7 9 2 8];
M_.block_structure_stat.equation_reordered = [ 9 10 3 4 5 6 7 8 1 2];
M_.block_structure_stat.incidence.sparse_IM = [
 1 2;
 1 8;
 1 10;
 2 3;
 2 4;
 2 8;
 2 10;
 3 1;
 3 4;
 3 5;
 3 7;
 4 3;
 4 4;
 4 7;
 5 2;
 5 5;
 5 7;
 6 5;
 6 6;
 7 6;
 7 7;
 7 8;
 7 9;
 8 8;
 8 9;
 8 10;
 9 1;
 10 10;
];
M_.block_structure_stat.tmp_nbr = 2;
M_.block_structure_stat.block(1).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(2).g1_sparse_rowval = int32([]);
M_.block_structure_stat.block(2).g1_sparse_colval = int32([]);
M_.block_structure_stat.block(2).g1_sparse_colptr = int32([]);
M_.block_structure_stat.block(3).g1_sparse_rowval = int32([1 2 8 2 8 1 3 4 4 5 1 2 3 5 5 6 3 7 5 6 7 8 ]);
M_.block_structure_stat.block(3).g1_sparse_colval = int32([1 1 1 2 2 3 3 3 4 4 5 5 5 5 6 6 7 7 8 8 8 8 ]);
M_.block_structure_stat.block(3).g1_sparse_colptr = int32([1 4 6 9 11 15 17 19 23 ]);
M_.static_g1_sparse_rowval = int32([3 9 1 5 2 4 2 3 4 3 5 6 6 7 3 4 5 7 1 2 7 8 7 8 1 2 8 10 ]);
M_.static_g1_sparse_colval = int32([1 1 2 2 3 3 4 4 4 5 5 5 6 6 7 7 7 7 8 8 8 8 9 9 10 10 10 10 ]);
M_.static_g1_sparse_colptr = int32([1 3 5 7 10 13 15 19 23 25 29 ]);
load  dynare_inputs; 
set_param_value('Phi'  ,dynare_inputs.param_values.Phi  );
set_param_value('delta' ,dynare_inputs.param_values.delta );
set_param_value('sigma' ,dynare_inputs.param_values.sigma );
set_param_value('alpha'   ,dynare_inputs.param_values.alpha   ); 
set_param_value('rho'   ,dynare_inputs.param_values.rho   );
set_param_value('chi'   ,dynare_inputs.param_values.chi   ); 
set_param_value('Z_ss'  ,dynare_inputs.ss(1)  );
set_param_value('R_ss'  ,dynare_inputs.ss(2) );
set_param_value('W_ss'  ,dynare_inputs.ss(3) );
set_param_value('L_ss'  ,dynare_inputs.ss(4)   );
set_param_value('K_ss'  ,dynare_inputs.ss(5)   ); 
set_param_value('I_ss'  ,dynare_inputs.ss(6)   );
set_param_value('Y_ss'  ,dynare_inputs.ss(7)   );
set_param_value('C_ss'  ,dynare_inputs.ss(8)   );
set_param_value('G_ss'  ,dynare_inputs.ss(9)   );
set_param_value('T_ss'  ,dynare_inputs.ss(10)   );
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(1) = M_.params(7);
oo_.steady_state(7) = M_.params(13);
oo_.steady_state(5) = M_.params(11);
oo_.steady_state(8) = M_.params(14);
oo_.steady_state(6) = M_.params(12);
oo_.steady_state(4) = M_.params(10);
oo_.steady_state(10) = M_.params(16);
oo_.steady_state(9) = M_.params(15);
oo_.steady_state(2) = M_.params(8);
oo_.steady_state(3) = M_.params(9);
if M_.exo_nbr > 0
	oo_.exo_simul = ones(M_.maximum_lag,1)*oo_.exo_steady_state';
end
if M_.exo_det_nbr > 0
	oo_.exo_det_simul = ones(M_.maximum_lag,1)*oo_.exo_det_steady_state';
end
%
% SHOCKS instructions
%
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = 0.1;
options_.impulse_responses.plot_threshold = 0;
options_.irf = 20;
options_.order = 1;
var_list_ = {'Z';'R';'W';'L';'K';'I';'Y';'C';'G';'T'};
[info, oo_, options_, M_] = stoch_simul(M_, options_, oo_, var_list_);


oo_.time = toc(tic0);
disp(['Total computing time : ' dynsec2hms(oo_.time) ]);
if ~exist([M_.dname filesep 'Output'],'dir')
    mkdir(M_.dname,'Output');
end
save([M_.dname filesep 'Output' filesep 'HPJedit_baseline_RBC_results.mat'], 'oo_', 'M_', 'options_');
if exist('estim_params_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'HPJedit_baseline_RBC_results.mat'], 'estim_params_', '-append');
end
if exist('bayestopt_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'HPJedit_baseline_RBC_results.mat'], 'bayestopt_', '-append');
end
if exist('dataset_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'HPJedit_baseline_RBC_results.mat'], 'dataset_', '-append');
end
if exist('estimation_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'HPJedit_baseline_RBC_results.mat'], 'estimation_info', '-append');
end
if exist('dataset_info', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'HPJedit_baseline_RBC_results.mat'], 'dataset_info', '-append');
end
if exist('oo_recursive_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'HPJedit_baseline_RBC_results.mat'], 'oo_recursive_', '-append');
end
if exist('options_mom_', 'var') == 1
  save([M_.dname filesep 'Output' filesep 'HPJedit_baseline_RBC_results.mat'], 'options_mom_', '-append');
end
if ~isempty(lastwarn)
  disp('Note: warning(s) encountered in MATLAB/Octave code')
end
