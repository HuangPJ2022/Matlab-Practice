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
M_.exo_names(1) = {'eps'};
M_.exo_names_tex(1) = {'eps'};
M_.exo_names_long(1) = {'eps'};
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
M_.endo_names(9) = {'V'};
M_.endo_names_tex(9) = {'V'};
M_.endo_names_long(9) = {'V'};
M_.endo_names(10) = {'d'};
M_.endo_names_tex(10) = {'d'};
M_.endo_names_long(10) = {'d'};
M_.endo_partitions = struct();
M_.param_names = cell(17,1);
M_.param_names_tex = cell(17,1);
M_.param_names_long = cell(17,1);
M_.param_names(1) = {'beta'};
M_.param_names_tex(1) = {'beta'};
M_.param_names_long(1) = {'beta'};
M_.param_names(2) = {'delta0'};
M_.param_names_tex(2) = {'delta0'};
M_.param_names_long(2) = {'delta0'};
M_.param_names(3) = {'sigma'};
M_.param_names_tex(3) = {'sigma'};
M_.param_names_long(3) = {'sigma'};
M_.param_names(4) = {'alpha'};
M_.param_names_tex(4) = {'alpha'};
M_.param_names_long(4) = {'alpha'};
M_.param_names(5) = {'rho'};
M_.param_names_tex(5) = {'rho'};
M_.param_names_long(5) = {'rho'};
M_.param_names(6) = {'omega'};
M_.param_names_tex(6) = {'omega'};
M_.param_names_long(6) = {'omega'};
M_.param_names(7) = {'gamma'};
M_.param_names_tex(7) = {'gamma'};
M_.param_names_long(7) = {'gamma'};
M_.param_names(8) = {'Z_ss'};
M_.param_names_tex(8) = {'Z\_ss'};
M_.param_names_long(8) = {'Z_ss'};
M_.param_names(9) = {'R_ss'};
M_.param_names_tex(9) = {'R\_ss'};
M_.param_names_long(9) = {'R_ss'};
M_.param_names(10) = {'W_ss'};
M_.param_names_tex(10) = {'W\_ss'};
M_.param_names_long(10) = {'W_ss'};
M_.param_names(11) = {'L_ss'};
M_.param_names_tex(11) = {'L\_ss'};
M_.param_names_long(11) = {'L_ss'};
M_.param_names(12) = {'K_ss'};
M_.param_names_tex(12) = {'K\_ss'};
M_.param_names_long(12) = {'K_ss'};
M_.param_names(13) = {'I_ss'};
M_.param_names_tex(13) = {'I\_ss'};
M_.param_names_long(13) = {'I_ss'};
M_.param_names(14) = {'Y_ss'};
M_.param_names_tex(14) = {'Y\_ss'};
M_.param_names_long(14) = {'Y_ss'};
M_.param_names(15) = {'C_ss'};
M_.param_names_tex(15) = {'C\_ss'};
M_.param_names_long(15) = {'C_ss'};
M_.param_names(16) = {'V_ss'};
M_.param_names_tex(16) = {'V\_ss'};
M_.param_names_long(16) = {'V_ss'};
M_.param_names(17) = {'d_ss'};
M_.param_names_tex(17) = {'d\_ss'};
M_.param_names_long(17) = {'d_ss'};
M_.param_partitions = struct();
M_.exo_det_nbr = 0;
M_.exo_nbr = 1;
M_.endo_nbr = 10;
M_.param_nbr = 17;
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
 0 11 15;
 0 12 16;]';
M_.nstatic = 4;
M_.nfwrd   = 4;
M_.npred   = 2;
M_.nboth   = 0;
M_.nsfwrd   = 4;
M_.nspred   = 2;
M_.ndynamic   = 6;
M_.dynamic_tmp_nbr = [6; 4; 0; 0; ];
M_.equations_tags = {
  1 , 'name' , '1' ;
  2 , 'name' , 'L' ;
  3 , 'name' , 'Y' ;
  4 , 'name' , '4' ;
  5 , 'name' , 'K' ;
  6 , 'name' , 'Z' ;
  7 , 'name' , 'R' ;
  8 , 'name' , '8' ;
  9 , 'name' , 'W' ;
  10 , 'name' , 'd' ;
};
M_.mapping.Z.eqidx = [3 6 7 9 ];
M_.mapping.R.eqidx = [1 7 8 ];
M_.mapping.W.eqidx = [2 9 ];
M_.mapping.L.eqidx = [2 3 7 9 ];
M_.mapping.K.eqidx = [3 5 7 9 ];
M_.mapping.I.eqidx = [4 5 ];
M_.mapping.Y.eqidx = [3 4 ];
M_.mapping.C.eqidx = [1 2 4 ];
M_.mapping.V.eqidx = [1 3 7 8 9 10 ];
M_.mapping.d.eqidx = [1 5 10 ];
M_.mapping.eps.eqidx = [6 ];
M_.static_and_dynamic_models_differ = false;
M_.has_external_function = false;
M_.block_structure.time_recursive = false;
M_.block_structure.block(1).Simulation_Type = 1;
M_.block_structure.block(1).endo_nbr = 1;
M_.block_structure.block(1).mfs = 1;
M_.block_structure.block(1).equation = [ 6];
M_.block_structure.block(1).variable = [ 1];
M_.block_structure.block(1).is_linear = true;
M_.block_structure.block(1).NNZDerivatives = 2;
M_.block_structure.block(1).bytecode_jacob_cols_to_sparse = [1 2 ];
M_.block_structure.block(2).Simulation_Type = 8;
M_.block_structure.block(2).endo_nbr = 9;
M_.block_structure.block(2).mfs = 9;
M_.block_structure.block(2).equation = [ 2 3 4 9 5 1 7 8 10];
M_.block_structure.block(2).variable = [ 3 7 6 4 5 8 9 2 10];
M_.block_structure.block(2).is_linear = false;
M_.block_structure.block(2).NNZDerivatives = 31;
M_.block_structure.block(2).bytecode_jacob_cols_to_sparse = [5 10 11 12 13 14 15 16 17 18 24 25 26 27 ];
M_.block_structure.block(1).g1_sparse_rowval = int32([]);
M_.block_structure.block(1).g1_sparse_colval = int32([]);
M_.block_structure.block(1).g1_sparse_colptr = int32([]);
M_.block_structure.block(2).g1_sparse_rowval = int32([2 4 5 7 1 4 2 3 3 5 1 2 4 7 5 1 3 6 2 4 7 8 9 7 8 5 9 6 6 6 6 ]);
M_.block_structure.block(2).g1_sparse_colval = int32([5 5 5 5 10 10 11 11 12 12 13 13 13 13 14 15 15 15 16 16 16 16 16 17 17 18 18 24 25 26 27 ]);
M_.block_structure.block(2).g1_sparse_colptr = int32([1 1 1 1 1 5 5 5 5 5 7 9 11 15 16 19 24 26 28 28 28 28 28 28 29 30 31 32 ]);
M_.block_structure.variable_reordered = [ 1 3 7 6 4 5 8 9 2 10];
M_.block_structure.equation_reordered = [ 6 2 3 4 9 5 1 7 8 10];
M_.block_structure.incidence(1).lead_lag = -1;
M_.block_structure.incidence(1).sparse_IM = [
 3 5;
 5 5;
 6 1;
 7 5;
 9 5;
];
M_.block_structure.incidence(2).lead_lag = 0;
M_.block_structure.incidence(2).sparse_IM = [
 1 8;
 2 3;
 2 4;
 2 8;
 3 1;
 3 4;
 3 7;
 3 9;
 4 6;
 4 7;
 4 8;
 5 5;
 5 6;
 5 10;
 6 1;
 7 1;
 7 2;
 7 4;
 7 9;
 8 2;
 8 9;
 9 1;
 9 3;
 9 4;
 9 9;
 10 9;
 10 10;
];
M_.block_structure.incidence(3).lead_lag = 1;
M_.block_structure.incidence(3).sparse_IM = [
 1 2;
 1 8;
 1 9;
 1 10;
];
M_.block_structure.dyn_tmp_nbr = 9;
M_.state_var = [1 5 ];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(10, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(1, 1);
M_.params = NaN(17, 1);
M_.endo_trends = struct('deflator', cell(10, 1), 'log_deflator', cell(10, 1), 'growth_factor', cell(10, 1), 'log_growth_factor', cell(10, 1));
M_.NNZDerivatives = [37; -1; -1; ];
M_.dynamic_g1_sparse_rowval = int32([6 3 5 7 9 3 6 7 9 7 8 2 9 2 3 7 9 5 4 5 3 4 1 2 4 3 7 8 9 10 5 10 1 1 1 1 6 ]);
M_.dynamic_g1_sparse_colval = int32([1 5 5 5 5 11 11 11 11 12 12 13 13 14 14 14 14 15 16 16 17 17 18 18 18 19 19 19 19 19 20 20 22 28 29 30 31 ]);
M_.dynamic_g1_sparse_colptr = int32([1 2 2 2 2 6 6 6 6 6 6 10 12 14 18 19 21 23 26 31 33 33 34 34 34 34 34 34 35 36 37 38 ]);
M_.lhs = {
'1/C'; 
'L'; 
'Y'; 
'Y'; 
'K'; 
'Z'; 
'R'; 
'R'; 
'W'; 
'd'; 
};
M_.static_tmp_nbr = [5; 4; 0; 0; ];
M_.block_structure_stat.block(1).Simulation_Type = 3;
M_.block_structure_stat.block(1).endo_nbr = 1;
M_.block_structure_stat.block(1).mfs = 1;
M_.block_structure_stat.block(1).equation = [ 6];
M_.block_structure_stat.block(1).variable = [ 1];
M_.block_structure_stat.block(2).Simulation_Type = 6;
M_.block_structure_stat.block(2).endo_nbr = 9;
M_.block_structure_stat.block(2).mfs = 9;
M_.block_structure_stat.block(2).equation = [ 2 3 4 5 1 7 8 9 10];
M_.block_structure_stat.block(2).variable = [ 3 5 7 6 8 9 2 4 10];
M_.block_structure_stat.variable_reordered = [ 1 3 5 7 6 8 9 2 4 10];
M_.block_structure_stat.equation_reordered = [ 6 2 3 4 5 1 7 8 9 10];
M_.block_structure_stat.incidence.sparse_IM = [
 1 2;
 1 8;
 1 9;
 1 10;
 2 3;
 2 4;
 2 8;
 3 1;
 3 4;
 3 5;
 3 7;
 3 9;
 4 6;
 4 7;
 4 8;
 5 5;
 5 6;
 5 10;
 6 1;
 7 1;
 7 2;
 7 4;
 7 5;
 7 9;
 8 2;
 8 9;
 9 1;
 9 3;
 9 4;
 9 5;
 9 9;
 10 9;
 10 10;
];
M_.block_structure_stat.tmp_nbr = 9;
M_.block_structure_stat.block(1).g1_sparse_rowval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colval = int32([1 ]);
M_.block_structure_stat.block(1).g1_sparse_colptr = int32([1 2 ]);
M_.block_structure_stat.block(2).g1_sparse_rowval = int32([1 8 2 4 6 8 2 3 3 4 1 3 5 2 5 6 7 8 9 5 6 7 1 2 6 8 4 5 9 ]);
M_.block_structure_stat.block(2).g1_sparse_colval = int32([1 1 2 2 2 2 3 3 4 4 5 5 5 6 6 6 6 6 6 7 7 7 8 8 8 8 9 9 9 ]);
M_.block_structure_stat.block(2).g1_sparse_colptr = int32([1 3 7 9 11 14 20 23 27 30 ]);
M_.static_g1_sparse_rowval = int32([3 6 7 9 1 7 8 2 9 2 3 7 9 3 5 7 9 4 5 3 4 1 2 4 1 3 7 8 9 10 1 5 10 ]);
M_.static_g1_sparse_colval = int32([1 1 1 1 2 2 2 3 3 4 4 4 4 5 5 5 5 6 6 7 7 8 8 8 9 9 9 9 9 9 10 10 10 ]);
M_.static_g1_sparse_colptr = int32([1 5 8 10 14 18 20 22 25 31 34 ]);
load  dynare_inputs; 
set_param_value('beta'  ,dynare_inputs.param_values.beta  );
set_param_value('delta0' ,dynare_inputs.param_values.delta0 );
set_param_value('sigma' ,dynare_inputs.param_values.sigma );
set_param_value('alpha'   ,dynare_inputs.param_values.alpha   ); 
set_param_value('rho'   ,dynare_inputs.param_values.rho   );
set_param_value('omega'   ,dynare_inputs.param_values.omega   ); 
set_param_value('gamma'   ,dynare_inputs.param_values.gamma   ); 
set_param_value('Z_ss'  ,dynare_inputs.ss_baseline(1)  );
set_param_value('R_ss'  ,dynare_inputs.ss_baseline(2) );
set_param_value('W_ss'  ,dynare_inputs.ss_baseline(3) );
set_param_value('L_ss'  ,dynare_inputs.ss_baseline(4)   );
set_param_value('K_ss'  ,dynare_inputs.ss_baseline(5)   ); 
set_param_value('I_ss'  ,dynare_inputs.ss_baseline(6)   );
set_param_value('Y_ss'  ,dynare_inputs.ss_baseline(7)   );
set_param_value('C_ss'  ,dynare_inputs.ss_baseline(8)   );
set_param_value('V_ss'  ,dynare_inputs.ss_baseline(9)   );
set_param_value('d_ss'  ,dynare_inputs.ss_baseline(10)   );
%
% INITVAL instructions
%
options_.initval_file = false;
oo_.steady_state(1) = M_.params(8);
oo_.steady_state(7) = M_.params(14);
oo_.steady_state(5) = M_.params(12);
oo_.steady_state(8) = M_.params(15);
oo_.steady_state(6) = M_.params(13);
oo_.steady_state(4) = M_.params(11);
oo_.steady_state(2) = M_.params(9);
oo_.steady_state(3) = M_.params(10);
oo_.steady_state(9) = M_.params(16);
oo_.steady_state(10) = M_.params(17);
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
M_.Sigma_e(1, 1) = 1;
options_.impulse_responses.plot_threshold = 0;
options_.irf = 20;
options_.order = 1;
var_list_ = {'Z';'R';'W';'L';'K';'I';'Y';'C';'V';'d'};
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
