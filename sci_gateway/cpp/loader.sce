// This file is released under the 3-clause BSD license. See COPYING-BSD.
// Generated by builder.sce : Please, do not edit this file
// ----------------------------------------------------------------------------
//
if ~win64() then
  warning(_("This module requires a Windows x64 platform."));
  return
end
//
FOSSEE_Optimization_path = get_absolute_file_path('loader.sce');
//
// ulink previous function with same name
[bOK, ilib] = c_link('FOSSEE_Optimization_Toolbox');
if bOK then
  ulink(ilib);
end
//
list_functions = [ 'sym_open';
                   'sym_close';
                   'sym_resetParams';
                   'sym_setIntParam';
                   'sym_getIntParam';
                   'sym_setDblParam';
                   'sym_getDblParam';
                   'sym_setStrParam';
                   'sym_getStrParam';
                   'sym_loadProblemBasic';
                   'sym_loadProblem';
                   'sym_solve';
                   'sym_getStatus';
                   'sym_isOptimal';
                   'sym_isInfeasible';
                   'sym_isAbandoned';
                   'sym_isIterLimitReached';
                   'sym_isTimeLimitReached';
                   'sym_isTargetGapAchieved';
                   'sym_getVarSoln';
                   'sym_getObjVal';
                   'sym_getIterCount';
                   'linearprog';
                   'rmps';
                   'solveqp';
                   'solveminuncp';
                   'solveminbndp';
                   'solveminconp';
                   'inter_fminunc';
                   'inter_fminbnd';
                   'inter_fmincon';
                   'sci_intqpipopt';
                   'sci_matrix_intlinprog';
                   'sci_mps_intlinprog';
                   'fotversion';
];
addinter(FOSSEE_Optimization_path + filesep() + 'FOSSEE_Optimization_Toolbox' + getdynlibext(), 'FOSSEE_Optimization_Toolbox', list_functions);
// remove temp. variables on stack
clear FOSSEE_Optimization_path;
clear bOK;
clear ilib;
clear list_functions;
// ----------------------------------------------------------------------------
