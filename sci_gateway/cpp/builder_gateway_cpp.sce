// Copyright (C) 2015 - IIT Bombay - FOSSEE
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt
// Author: Harpreet Singh
// Organization: FOSSEE, IIT Bombay
// Email: toolbox@scilab.in

mode(-1)
lines(0)

toolbox_title = "FOSSEE_Optimization_Toolbox";

Build_64Bits = %t;

path_builder = get_absolute_file_path('builder_gateway_cpp.sce');


if getos()=="Windows" then
//Name of All the Functions
Function_Names = [
		//for opening/closing environment and checking if it is open/close
		"sym_open","sci_sym_open";
		"sym_close","sci_sym_close";
		"sym_isEnvActive","sci_sym_isEnvActive";
		
		//run time parameters
		"sym_resetParams","sci_sym_set_defaults";
		"sym_setIntParam","sci_sym_set_int_param";
		"sym_getIntParam","sci_sym_get_int_param";
		"sym_setDblParam","sci_sym_set_dbl_param";
		"sym_getDblParam","sci_sym_get_dbl_param";
		"sym_setStrParam","sci_sym_set_str_param";
		"sym_getStrParam","sci_sym_get_str_param";
		"sym_getInfinity","sci_sym_getInfinity";
		
		//problem loaders
		"sym_loadProblemBasic","sci_sym_loadProblemBasic";
		"sym_loadProblem","sci_sym_loadProblem";
		"sym_loadMPS","sci_sym_load_mps";
		
		//basic data
		"sym_getNumConstr","sci_sym_get_num_int";
		"sym_getNumVar","sci_sym_get_num_int";
		"sym_getNumElements","sci_sym_get_num_int";
		
		//variable and objective data
		"sym_isContinuous","sci_sym_isContinuous";
		"sym_isBinary","sci_sym_isBinary";
		"sym_isInteger","sci_sym_isInteger";
		"sym_setContinuous","sci_sym_set_continuous";
		"sym_setInteger","sci_sym_set_integer";
		"sym_getVarLower","sci_sym_get_dbl_arr";
		"sym_getVarUpper","sci_sym_get_dbl_arr";
		"sym_setVarLower","sci_sym_setVarBound";
		"sym_setVarUpper","sci_sym_setVarBound";
		"sym_getObjCoeff","sci_sym_get_dbl_arr";
		"sym_setObjCoeff","sci_sym_setObjCoeff";
		"sym_getObjSense","sci_sym_getObjSense";
		"sym_setObjSense","sci_sym_setObjSense";
		
		//constraint data
		"sym_getRhs","sci_sym_get_dbl_arr";
		"sym_getConstrRange","sci_sym_get_dbl_arr";
		"sym_getConstrLower","sci_sym_get_dbl_arr";
		"sym_getConstrUpper","sci_sym_get_dbl_arr";
		"sym_setConstrLower","sci_sym_setConstrBound";
		"sym_setConstrUpper","sci_sym_setConstrBound";
		"sym_setConstrType","sci_sym_setConstrType";
		"sym_getMatrix","sci_sym_get_matrix";
		
		//add/remove variables and constraints
		"sym_addConstr","sci_sym_addConstr";
		"sym_addVar","sci_sym_addVar";
		"sym_deleteVars","sci_sym_delete_cols";
		"sym_deleteConstrs","sci_sym_delete_rows";
		
		//primal bound
		"sym_getPrimalBound","sci_sym_getPrimalBound";
		"sym_setPrimalBound","sci_sym_setPrimalBound";
		
		//set preliminary solution
		"sym_setVarSoln","sci_sym_setColSoln";
		
		//solve
		"sym_solve","sci_sym_solve";
		
		//post solve functions
		"sym_getStatus","sci_sym_get_status";
		"sym_isOptimal","sci_sym_get_solver_status";
		"sym_isInfeasible","sci_sym_get_solver_status";
		"sym_isAbandoned","sci_sym_get_solver_status";
		"sym_isIterLimitReached","sci_sym_get_solver_status";
		"sym_isTimeLimitReached","sci_sym_get_solver_status";
		"sym_isTargetGapAchieved","sci_sym_get_solver_status";
		"sym_getVarSoln","sci_sym_getVarSoln";
		"sym_getObjVal","sci_sym_getObjVal";
		"sym_getIterCount","sci_sym_get_iteration_count";
		"sym_getConstrActivity","sci_sym_getRowActivity";

		//Linprog function
		"linearprog","sci_linearprog"
        "rmps","sci_rmps"

		//QP function
		"solveqp","sci_solveqp"

		//fminunc function and fminbnd function
		"solveminuncp","sci_solveminuncp"
		"solveminbndp","sci_solveminbndp"
		"solveminconp","sci_solveminconp"
		
	];

//Name of all the files to be compiled
Files = [
		"globals.cpp",
		"sci_iofunc.cpp",
		"sci_sym_openclose.cpp",
		"sci_solver_status_query_functions.cpp",
		"sci_sym_solve.cpp",                    
		"sci_sym_loadproblem.cpp",
		"sci_sym_isenvactive.cpp",
		"sci_sym_load_mps.cpp",
		"sci_vartype.cpp",
		"sci_sym_getinfinity.cpp",
		"sci_sym_solution.cpp",
		"sci_sym_get_dbl_arr.cpp",
		"sci_sym_get_iteration_count.cpp",
		"sci_sym_get_matrix.cpp",
		"sci_sym_get_num_int.cpp",
		"sci_sym_set_variables.cpp",
		"sci_sym_setobj.cpp",
		"sci_sym_varbounds.cpp",
		"sci_sym_rowmod.cpp",
		"sci_sym_set_indices.cpp",
		"sci_sym_addrowcol.cpp",
		"sci_sym_primalbound.cpp",
		"sci_sym_setcolsoln.cpp",
		"sci_sym_getrowact.cpp",
		"sci_sym_getobjsense.cpp",
		"sci_sym_remove.cpp",
		"sci_QuadNLP.cpp",
		"sci_ipopt.cpp",
		"sci_QuadNLP.cpp",
		"sci_ipopt.cpp",
		"sci_minuncNLP.cpp",
		"sci_ipoptfminunc.cpp",
		"sci_minbndNLP.cpp",
		"sci_ipoptfminbnd.cpp",
		"sci_minconNLP.cpp",
		"sci_ipoptfmincon.cpp",
		"sci_LinProg.cpp",
        "read_mps.cpp"
	]
else
//Name of All the Functions
Function_Names = [
	

		//QP(IPOPT)
		"solveqp","sci_solveqp";

		//fmincon,fminunc, fminbnd function(IPOPT)
		"solveminuncp","sci_solveminuncp";
		"solveminbndp","sci_solveminbndp";
		"solveminconp","sci_solveminconp";
		
		//Linprog function (CLP)
		"linearprog","sci_linearprog";
		"rmps", "sci_rmps";

		//MILP (CBC)
		'sci_matrix_intlinprog', 'matrix_cppintlinprog';
		'sci_mps_intlinprog','mps_cppintlinprog';

				
		//Integer NLP (Bonmin)
		'inter_fminunc', 'cpp_intfminunc';
		'inter_fminbnd', 'cpp_intfminbnd';
		'inter_fmincon', 'cpp_intfmincon';
		'sci_intqpipopt', 'cpp_intqpipopt';

		//Symphony

		//for opening/closing environment and checking if it is open/close
		"sym_open","sci_sym_open";
		"sym_close","sci_sym_close";
		"sym_isEnvActive","sci_sym_isEnvActive";

		//run time parameters
		"sym_resetParams","sci_sym_set_defaults";
		"sym_setIntParam","sci_sym_set_int_param";
		"sym_getIntParam","sci_sym_get_int_param";
		"sym_setDblParam","sci_sym_set_dbl_param";
		"sym_getDblParam","sci_sym_get_dbl_param";
		"sym_setStrParam","sci_sym_set_str_param";
		"sym_getStrParam","sci_sym_get_str_param";
		"sym_getInfinity","sci_sym_getInfinity";

		//problem loaders
		"sym_loadProblemBasic","sci_sym_loadProblemBasic";
		"sym_loadProblem","sci_sym_loadProblem";
		"sym_loadMPS","sci_sym_load_mps";

		//basic data
		"sym_getNumConstr","sci_sym_get_num_int";
		"sym_getNumVar","sci_sym_get_num_int";
		"sym_getNumElements","sci_sym_get_num_int";

		//variable and objective data
		"sym_isContinuous","sci_sym_isContinuous";
		"sym_isBinary","sci_sym_isBinary";
		"sym_isInteger","sci_sym_isInteger";
		"sym_setContinuous","sci_sym_set_continuous";
		"sym_setInteger","sci_sym_set_integer";
		"sym_getVarLower","sci_sym_get_dbl_arr";
		"sym_getVarUpper","sci_sym_get_dbl_arr";
		"sym_setVarLower","sci_sym_setVarBound";
		"sym_setVarUpper","sci_sym_setVarBound";
		"sym_getObjCoeff","sci_sym_get_dbl_arr";
		"sym_setObjCoeff","sci_sym_setObjCoeff";
		"sym_getObjSense","sci_sym_getObjSense";
		"sym_setObjSense","sci_sym_setObjSense";

		//constraint data
		"sym_getRhs","sci_sym_get_dbl_arr";
		"sym_getConstrRange","sci_sym_get_dbl_arr";
		"sym_getConstrLower","sci_sym_get_dbl_arr";
		"sym_getConstrUpper","sci_sym_get_dbl_arr";
		"sym_setConstrLower","sci_sym_setConstrBound";
		"sym_setConstrUpper","sci_sym_setConstrBound";
		"sym_setConstrType","sci_sym_setConstrType";
		"sym_getMatrix","sci_sym_get_matrix";

		//add/remove variables and constraints
		"sym_addConstr","sci_sym_addConstr";
		"sym_addVar","sci_sym_addVar";
		"sym_deleteVars","sci_sym_delete_cols";
		"sym_deleteConstrs","sci_sym_delete_rows";

		//primal bound
		"sym_getPrimalBound","sci_sym_getPrimalBound";
		"sym_setPrimalBound","sci_sym_setPrimalBound";
		
		//set preliminary solution
		"sym_setVarSoln","sci_sym_setColSoln";
	
		//solve
		"sym_solve","sci_sym_solve";
		
		//post solve functions
		"sym_getStatus","sci_sym_get_status";
		"sym_isOptimal","sci_sym_get_solver_status";
		"sym_isInfeasible","sci_sym_get_solver_status";
		"sym_isAbandoned","sci_sym_get_solver_status";
		"sym_isIterLimitReached","sci_sym_get_solver_status";
		"sym_isTimeLimitReached","sci_sym_get_solver_status";
		"sym_isTargetGapAchieved","sci_sym_get_solver_status";
		"sym_getVarSoln","sci_sym_getVarSoln";
		"sym_getObjVal","sci_sym_getObjVal";
		"sym_getIterCount","sci_sym_get_iteration_count";
		"sym_getConstrActivity","sci_sym_getRowActivity";

		
		//ecos function
		"solveecos","sci_ecos";
   

		"fotversion","sci_fotversion"
	];

//Name of all the files to be compiled
Files = [
		
		"sci_iofunc.cpp",


		//QP (IPOPT)
		"sci_ipopt.cpp",
		"sci_QuadNLP.cpp",
	
		//fmincon,fminunc,fminbnd (IPOPT)
		"sci_minuncNLP.cpp",
		"sci_ipoptfminunc.cpp",
		"sci_minbndNLP.cpp",
		"sci_ipoptfminbnd.cpp",
		"sci_minconNLP.cpp",
		"sci_ipoptfmincon.cpp",


		//Linprog (CLP)
		"sci_LinProg.cpp",
	     "read_mps.cpp",

		//MILP (CBC)
		'cbcintlinprog_mpscpp.cpp',
		'cbcintlinprog_matrixcpp.cpp',
		
		//Integer NLP (Bonmin) 
		'sci_QuadTMINLP.cpp',
		'cpp_intqpipopt.cpp',
		'sci_minuncTMINLP.cpp',
		'cpp_intfminunc.cpp',
		'sci_minbndTMINLP.cpp',
		'cpp_intfminbnd.cpp',		
		'sci_minconTMINLP.cpp',
		'cpp_intfmincon.cpp',

	
		//Symphony
		"globals.cpp",
		"sci_sym_openclose.cpp",
		"sci_solver_status_query_functions.cpp",
		"sci_sym_solve.cpp",                    
		"sci_sym_loadproblem.cpp",
		"sci_sym_isenvactive.cpp",
		"sci_sym_load_mps.cpp",
		"sci_vartype.cpp",
		"sci_sym_getinfinity.cpp",
		"sci_sym_solution.cpp",
		"sci_sym_get_dbl_arr.cpp",
		"sci_sym_get_iteration_count.cpp",
		"sci_sym_get_matrix.cpp",
		"sci_sym_get_num_int.cpp",
		"sci_sym_set_variables.cpp",
		"sci_sym_setobj.cpp",
		"sci_sym_varbounds.cpp",
		"sci_sym_rowmod.cpp",
		"sci_sym_set_indices.cpp",
		"sci_sym_addrowcol.cpp",
		"sci_sym_primalbound.cpp",
		"sci_sym_setcolsoln.cpp",
		"sci_sym_getrowact.cpp",
		"sci_sym_getobjsense.cpp",
		"sci_sym_remove.cpp", 

		//ecos
 		"ecos.cpp",

		"sci_fotversion.cpp"
	]

end


[a, opt] = getversion();
Version = opt(2);

//Build_64Bits = %f;

if getos()=="Windows" then
	third_dir = path_builder+filesep()+'..'+filesep()+'..'+filesep()+'thirdparty';
	lib_base_dir = third_dir + filesep() + 'windows' + filesep() + 'lib' + filesep() + Version + filesep();
	inc_base_dir = third_dir + filesep() + 'windows' + filesep() + 'include' + filesep() + 'coin';
    C_Flags=['-D__USE_DEPRECATED_STACK_FUNCTIONS__ -w -I '+path_builder+' '+ '-I '+inc_base_dir+' ']
    Linker_Flag  = [lib_base_dir+"libClp.lib "+lib_base_dir+"libCgl.lib "+lib_base_dir+"libOsi.lib "+lib_base_dir+"libOsiClp.lib "+lib_base_dir+"libCoinUtils.lib "+lib_base_dir+"libSymphony.lib "+lib_base_dir+"IpOptFSS.lib "+lib_base_dir+"IpOpt-vc10.lib "]

else
	third_dir = path_builder+filesep()+'..'+filesep()+'..'+filesep()+'thirdparty';
	lib_base_dir = third_dir + filesep() + 'linux' + filesep() + 'lib' + filesep() + Version + filesep();
	inc_base_dir = third_dir + filesep() + 'linux' + filesep() + 'include' + filesep() + 'coin'+ filesep();
    
    C_Flags=["-g -D__USE_DEPRECATED_STACK_FUNCTIONS__ -w -fpermissive -I"+path_builder+" -I"+inc_base_dir+" -Wl,-rpath="+lib_base_dir+" "]
    
    Linker_Flag = ["-L"+lib_base_dir+"libSym"+" "+"-L"+lib_base_dir+"libipopt"+" "+"-L"+lib_base_dir+"libClp"+" "+"-L"+lib_base_dir+"libOsiClp"+" "+"-L"+lib_base_dir+"libCoinUtils" ]

    
end

tbx_build_gateway(toolbox_title,Function_Names,Files,path_builder, [], Linker_Flag, C_Flags);


clear toolbox_title Function_Names Files Linker_Flag C_Flags;
