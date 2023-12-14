#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include "util_String.h"
#include "util_ErrorCodes.h"
#include "util_Table.h"

#include "EOS_Utils.h"

#define EOSBASE_DBG 1

/*
  This part is for thorns providing an EOS.
*/

/*
  The structure for a single EOS function.
*/

typedef struct
{
  const char* eos_name; /* The name of the EOS, required to match the call */
  CCTK_INT param_table; /* The parameter table allows the EOS 
                           to store persistent information */
  CCTK_INT n_indeps;    /* The number of independent variables 
                           the EOS requires */
  CCTK_INT n_deps;      /* The number of   dependent variables 
                           the EOS can return */
  char **indep_names;   /* The names (labels) of the independent variables */
  char **dep_names;     /* The names (labels) of the   dependent variables */
  CCTK_INT (*EOS_fn)(const CCTK_INT, 
                     const CCTK_INT,
                     const CCTK_REAL**,
                     const CCTK_INT*,
                     CCTK_REAL**); /* The actual function pointer */
} EOS_Function;

/*
  The list of all EOS functions.
*/

struct EOSFunctionList
{
  CCTK_INT n_eos_fns;
  EOS_Function* Fns;
} FnList;

 /*@@
   @routine    SetupFnList
   @date       Fri Jun 10 04:37:10 2005
   @author     Ian Hawke
   @desc 
   Initialize the list of functions.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void SetupFnList()
{
  FnList.n_eos_fns = 0;
  FnList.Fns = NULL;
}

 /*@@
   @routine    EOSBase_Register
   @date       Fri Jun 10 04:37:51 2005
   @author     Ian Hawke
   @desc 
   The (aliased) function called by an EOS thorn to register its call.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

CCTK_INT 
EOSBase_Register(const CCTK_INT registered_table_handle,
                 CCTK_INT registered_EOS_fn(const CCTK_INT param_table, 
                                            const CCTK_INT n_dims,
                                            const CCTK_POINTER* indep_vars,
                                            const CCTK_INT* which_deps_to_set,
                                            CCTK_POINTER* dep_vars))
{

  CCTK_INT ierr= 0, nfn, var, table_ierr;
  const char *tmp_indep_names;
  const char *tmp_dep_names;
  const char* tmpbuffer;

  CCTK_INT registered_n_indeps, registered_n_deps, strlength;

  ++(FnList.n_eos_fns);
  FnList.Fns = (EOS_Function*)realloc(FnList.Fns,
                                      (FnList.n_eos_fns) * 
                                      sizeof(EOS_Function));
  nfn = FnList.n_eos_fns - 1;

  (FnList.Fns)[nfn].param_table = registered_table_handle;

  table_ierr = Util_TableGetInt(registered_table_handle,
                                &registered_n_indeps,
                                "N independent variables");
  if (table_ierr != 1)
  {
    CCTK_WARN(0, "What's the number of independent variables?");
  }
  
  table_ierr = Util_TableGetInt(registered_table_handle,
                                &registered_n_deps,
                                "N dependent variables");
  if (table_ierr != 1)
  {
    CCTK_WARN(0, "What's the number of dependent variables?");
  }

  (FnList.Fns)[nfn].n_indeps = registered_n_indeps;
  (FnList.Fns)[nfn].n_deps = registered_n_deps;

  strlength = Util_TableGetString(registered_table_handle,
                                  0, NULL,
                                  "EOS Name");
  if (strlength < 1)
  {
    CCTK_WARN(0, "What EOS do you actually want?");
  }
  (FnList.Fns)[nfn].eos_name = 
    (const char*)malloc((strlength+1) * sizeof(char));
  ierr = Util_TableGetString(registered_table_handle,
                             (strlength+1), (FnList.Fns)[nfn].eos_name,
                             "EOS name");
  if (ierr != strlength)
  {
    CCTK_WARN(0, "What happened here?");
  }

  strlength = Util_TableGetString(registered_table_handle,
                                  0, NULL,
                                  "Independent variable names");
  if (strlength < 1)
  {
    CCTK_WARN(0, "Which independent do you actually want?");
  }
  tmp_indep_names = (const char*)malloc((strlength+1) * sizeof(char));
  ierr = Util_TableGetString(registered_table_handle,
                             (strlength+1), tmp_indep_names,
                             "Independent variable names");
  if (ierr != strlength)
  {
    CCTK_WARN(0, "What happened here?");
  }

  strlength = Util_TableGetString(registered_table_handle,
                                  0, NULL,
                                  "Dependent variable names");
  if (strlength < 1)
  {
    CCTK_WARN(0, "Which dependent variables do you actually want?");
  }
  tmp_dep_names = (const char*)malloc((strlength+1) * sizeof(char));
  ierr = Util_TableGetString(registered_table_handle,
                             (strlength+1), tmp_dep_names,
                             "Dependent variable names");
  if (ierr != strlength)
  {
    CCTK_WARN(0, "What happened here?");
  }
  
  (FnList.Fns)[nfn].EOS_fn = registered_EOS_fn;

  (FnList.Fns)[nfn].indep_names = (char**)malloc(registered_n_indeps * 
                                                 sizeof(char*));
  tmpbuffer = EOS_StrSep(&tmp_indep_names, " ");
  if (tmpbuffer == NULL) /* line is empty */      
  {
    CCTK_WARN(0, "Error in the independent name string; empty???");
  }
  for (var = 0; var < registered_n_indeps; ++var)
  {
    while (*tmpbuffer == '\0')
    {
      tmpbuffer = EOS_StrSep(&tmp_indep_names, " ");
    }
    ((FnList.Fns)[nfn].indep_names)[var] = Util_Strdup(tmpbuffer);
    tmpbuffer = EOS_StrSep(&tmp_indep_names, " ");
  }
  
  (FnList.Fns)[nfn].dep_names = (char**)malloc(registered_n_deps * 
                                               sizeof(char*));;
  tmpbuffer = EOS_StrSep(&tmp_dep_names, " ");
  if (tmpbuffer == NULL) /* line is empty */      
  {
    CCTK_WARN(0, "Error in the dependent name string; empty???");
  }
  for (var = 0; var < registered_n_deps; ++var)
  {
    while (*tmpbuffer == '\0')
    {
      tmpbuffer = EOS_StrSep(&tmp_dep_names, " ");
    }
    ((FnList.Fns)[nfn].dep_names)[var] = Util_Strdup(tmpbuffer);
    tmpbuffer = EOS_StrSep(&tmp_dep_names, " ");
  }
  
  ierr = 0;

#ifdef EOSBASE_DBG
  printf("An eos call has been registered. The call number was %d.\n"
         "The eos name is '%s'.\n"
         "The number of independent variables is %d.\n"
         "The number of dependent variables is %d.\n",
         FnList.n_eos_fns,
         (FnList.Fns)[nfn].eos_name,
         (FnList.Fns)[nfn].n_indeps,
         (FnList.Fns)[nfn].n_deps);
  for (var = 0; var < (FnList.Fns)[nfn].n_indeps; ++var)
  {
    printf("Independent variable %d is '%s'\n",
           var, (FnList.Fns)[nfn].indep_names[var]);
  }
  for (var = 0; var < (FnList.Fns)[nfn].n_deps; ++var)
  {
    printf("Dependent variable %d is '%s'\n",
           var, (FnList.Fns)[nfn].dep_names[var]);
  }
#endif
  
  return ierr;
}

/*
  This part is for thorns that want to call the EOS.
*/

/*
  The structure describing a single EOS 'call'.
*/

typedef struct
{
  CCTK_INT fn_number;  /* The number (in the EOS_FnList) of the EOS function */
  CCTK_INT* indep_gfs; /* List of independent GF indices
                          (if not using arrays) */
  CCTK_INT* dep_gfs;   /* List of   dependent GF indices
                          (if not using arrays) */
  CCTK_INT* which_deps_to_set; /* List of which of the dependent variables
                                  should be set */
  CCTK_INT n_elems;       /* Number of elements in the GFs/arrays */
  CCTK_REAL** indep_vars; /* Pointers to the independent variable arrays 
                             (if not using GFs) */
  CCTK_REAL** dep_vars;   /* Pointers to the   dependent variable arrays 
                             (if not using GFs) */
} EOS_Call;

/*
  The structure describing the list of all EOS calls.
*/

struct EOSCallList
{
  CCTK_INT n_eos_calls;
  EOS_Call* Calls;
} CallList;

 /*@@
   @routine    SetupCallList
   @date       Fri Jun 10 04:46:17 2005
   @author     Ian Hawke
   @desc 
   Initialize the call list.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void SetupCallList()
{
  CallList.n_eos_calls = 0;
  CallList.Calls = NULL;
}

 /*@@
   @routine    SetupCall
   @date       Fri Jun 10 04:46:31 2005
   @author     Ian Hawke
   @desc 
   Set up a single function call.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

CCTK_INT SetupCall(CCTK_INT table_handle)
{

  CCTK_INT ncall, i, strlength, call_number, type_code, n_entries, 
    ierr, n_elements, var, indep_var, dep_var;
  const char* fn_name;
  const char *local_indep_names;
  const char *local_dep_names;
  char *tmpbuffer;
  CCTK_INT *tmp_indep_gfs;
  CCTK_INT *tmp_dep_gfs;
  CCTK_REAL **tmp_indep_arrays;
  CCTK_REAL **tmp_dep_arrays;

  ++(CallList.n_eos_calls);
  CallList.Calls = (EOS_Call*)realloc(CallList.Calls,
                                      (CallList.n_eos_calls) * 
                                      sizeof(EOS_Call));
  ncall = CallList.n_eos_calls - 1;
  
  strlength = Util_TableGetString(table_handle,
                                  0, NULL,
                                  "EOS Name");
  if (strlength < 1)
  {
    CCTK_WARN(0, "What EOS do you actually want?");
  }
  fn_name = (const char*)malloc((strlength+1) * sizeof(char));
  ierr = Util_TableGetString(table_handle,
                             (strlength+1), fn_name,
                             "EOS Name");
  if (ierr != strlength)
  {
    CCTK_WARN(0, "What happened here?");
  }
  
  for (call_number = 0; call_number < FnList.n_eos_fns; ++call_number)
  {
    if (CCTK_Equals(fn_name, FnList.Fns[call_number].eos_name)) break;
  }
  if (call_number < FnList.n_eos_fns)
  {
    CallList.Calls[ncall].fn_number = call_number;
  }
  else
  {
    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "Failed to find eos function '%s'.",
               fn_name);
  }

  if (Util_TableQueryValueInfo(table_handle, &type_code, &n_entries,
                               "Independent GFs"))
  {
    if (type_code != CCTK_VARIABLE_INT)
    {
      CCTK_WARN(0, "The independent GF indices must be passed as "
                "a CCTK_INT array.");
    }
    if (n_entries != (FnList.Fns[call_number]).n_indeps)
    {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "The EOS function '%s' expects %d independent variables.\n"
                 "Your GF list contains %d entries.",
                 fn_name, (FnList.Fns[call_number]).n_indeps,
                 n_entries);
    }
    tmp_indep_gfs = (CCTK_INT*)malloc(n_entries * sizeof(CCTK_INT));
    (CallList.Calls[ncall]).indep_gfs = (CCTK_INT*)malloc(n_entries *
                                                          sizeof(CCTK_INT));
    n_elements = Util_TableGetIntArray(table_handle, 
                                       (FnList.Fns[call_number]).n_indeps,
                                       tmp_indep_gfs,
                                       "Independent GFs");
    if (n_elements != n_entries)
    {
      CCTK_WARN(0, "What happened here?");
    }
    /* Now we have to order the independent gfs correctly */
    
    strlength = Util_TableGetString(table_handle,
                                    0, NULL,
                                    "Independent variable names");
    if (strlength < 1)
    {
      CCTK_WARN(0, "Need the independent variable names");
    }
    local_indep_names = 
      (const char*)malloc((strlength+1) * sizeof(char));
    ierr = Util_TableGetString(table_handle,
                               (strlength+1), local_indep_names,
                               "Independent variable names");
    if (ierr != strlength)
    {
      CCTK_WARN(0, "What happened here?");
    }

    tmpbuffer = EOS_StrSep(&local_indep_names, " ");
    if (tmpbuffer == NULL) /* line is empty */      
    {
      CCTK_WARN(0, "Error in the independent name string; empty???");
    }
    for (var = 0; var < n_entries; ++var)
    {
      while (*tmpbuffer == '\0')
      {
        tmpbuffer = EOS_StrSep(&local_indep_names, " ");
      }

      (CallList.Calls[ncall]).indep_gfs[var] = -1;
      for (indep_var = 0; indep_var < n_entries; ++indep_var)
      {
        if ( CCTK_Equals((FnList.Fns[call_number]).indep_names[indep_var],
                         tmpbuffer) )
        {
          (CallList.Calls[ncall]).indep_gfs[var] = tmp_indep_gfs[indep_var];
        }
      }
      if ( (CallList.Calls[ncall]).indep_gfs[var] < 0 )
      {
        CCTK_WARN(0, "Failed to set up the independent variables?");
      }
      tmpbuffer = EOS_StrSep(&local_indep_names, " ");
    }

    if (Util_TableQueryValueInfo(table_handle, &type_code, &n_entries,
                                 "Dependent GFs"))
    {
      if (type_code != CCTK_VARIABLE_INT)
      {
        CCTK_WARN(0, "The dependent GF indices must be passed as "
                  "a CCTK_INT array.");
      }
      if (n_entries > (FnList.Fns[call_number]).n_deps)
      {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "The EOS function '%s' expects at most %d dependent"
                   " variables.\n"
                   "Your GF list contains %d entries.",
                   fn_name, (FnList.Fns[call_number]).n_deps,
                   n_entries);
      }
      tmp_dep_gfs = (CCTK_INT*)malloc(n_entries * sizeof(CCTK_INT));
      (CallList.Calls[ncall]).dep_gfs = 
        (CCTK_INT*)malloc((FnList.Fns[call_number]).n_deps * sizeof(CCTK_INT));
      (CallList.Calls[ncall]).which_deps_to_set = 
        (CCTK_INT*)malloc((FnList.Fns[call_number]).n_deps * sizeof(CCTK_INT));
      n_elements = Util_TableGetIntArray(table_handle, 
                                         n_entries,
                                         tmp_dep_gfs,
                                         "Dependent GFs");
      if (n_elements != n_entries)
      {
        CCTK_WARN(0, "What happened here?");
      }
      
      /* Now we have to order the dependent gfs correctly */
    
      strlength = Util_TableGetString(table_handle,
                                      0, NULL,
                                      "Dependent variable names");
      if (strlength < 1)
      {
        CCTK_WARN(0, "Need the dependent variable names");
      }
      local_dep_names = 
        (const char*)malloc((strlength+1) * sizeof(char));
      ierr = Util_TableGetString(table_handle,
                                 (strlength+1), local_dep_names,
                                 "Dependent variable names");
      if (ierr != strlength)
      {
        CCTK_WARN(0, "What happened here?");
      }

      tmpbuffer = EOS_StrSep(&local_dep_names, " ");
      if (tmpbuffer == NULL) /* line is empty */      
      {
        CCTK_WARN(0, "Error in the independent name string; empty???");
      }
      for (dep_var = 0; 
           dep_var < (FnList.Fns[call_number]).n_deps; 
           ++dep_var)
      {
        (CallList.Calls[ncall]).dep_gfs[dep_var] = -1;
        (CallList.Calls[ncall]).which_deps_to_set[dep_var] = 0;
      }
      for (var = 0; var < n_entries; ++var)
      {
        while (*tmpbuffer == '\0')
        {
          tmpbuffer = EOS_StrSep(&local_dep_names, " ");
        }

        for (dep_var = 0; 
             dep_var < (FnList.Fns[call_number]).n_deps; 
             ++dep_var)
        {
          if ( CCTK_Equals((FnList.Fns[call_number]).dep_names[dep_var],
                           tmpbuffer) )
          {
            (CallList.Calls[ncall]).dep_gfs[dep_var] = tmp_dep_gfs[var];
            (CallList.Calls[ncall]).which_deps_to_set[dep_var] = 1;
          }
        }
        tmpbuffer = EOS_StrSep(&local_dep_names, " ");
      }
    }
    else
    {
      CCTK_WARN(0, "If you set the 'Independent GFs' key you should also\n"
                "set the 'Dependent GFs' key.");
    }
    (CallList.Calls[ncall]).indep_vars = 
      (CCTK_REAL**)malloc((FnList.Fns[call_number]).n_indeps *
                                       sizeof(CCTK_REAL*));
    for (var = 0; var < (FnList.Fns[call_number]).n_indeps; ++var)
    {
      (CallList.Calls[ncall]).indep_vars[var] = NULL;
    }
    (CallList.Calls[ncall]).dep_vars = 
      (CCTK_REAL**)malloc((FnList.Fns[call_number]).n_deps *
                                     sizeof(CCTK_REAL*));
    for (var = 0; var < (FnList.Fns[call_number]).n_deps; ++var)
    {
      (CallList.Calls[ncall]).dep_vars[var] = NULL;
    }
  }
  else /* Setting up a direct array call */
  {
    (CallList.Calls[ncall]).indep_gfs = NULL;
    (CallList.Calls[ncall]).dep_gfs = NULL;
   
    if (Util_TableQueryValueInfo(table_handle, &type_code, &n_entries,
                                 "N array elements"))
    {
      if (type_code != CCTK_VARIABLE_INT)
      {
        CCTK_WARN(0, "The number of elements of the array in a direct "
                  "array call must be passed as a CCTK_INT.");
      }
      n_elements = Util_TableGetInt(table_handle,
                                    &((CallList.Calls[ncall]).n_elems),
                                    "N array elements");
    }
    else
    {
      CCTK_WARN(0, "If setting up a direct array call the number of "
                "elements\n(key 'N array elements') must be set");
    }
    
    if (Util_TableQueryValueInfo(table_handle, &type_code, &n_entries,
                                 "Independent arrays"))
    {
      if (type_code != CCTK_VARIABLE_POINTER)
      {
        CCTK_WARN(0, "The independent arrays in a direct "
                  "array call must be passed as as a CCTK_POINTER array.");
      }
      if (n_entries != (FnList.Fns[call_number]).n_indeps)
      {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "The EOS function '%s' expects %d independent variables.\n"
                   "Your direct call list contains %d entries.",
                   fn_name, (FnList.Fns[call_number]).n_indeps,
                   n_entries);
      }
      tmp_indep_arrays = (CCTK_REAL**)malloc(n_entries * sizeof(CCTK_REAL*));
      (CallList.Calls[ncall]).indep_vars = 
        (CCTK_REAL**)malloc(n_entries * sizeof(CCTK_REAL*));
      n_elements = 
        Util_TableGetPointerArray(table_handle,
                                  (FnList.Fns[call_number]).n_indeps,
                                  (CCTK_POINTER*)tmp_indep_arrays,
                                  "Independent arrays");
      if (n_elements != n_entries)
      {
        CCTK_WARN(0, "What happened here?");
      }
      /* Now we have to order the independent gfs correctly */
    
      strlength = Util_TableGetString(table_handle,
                                      0, NULL,
                                      "Independent variable names");
      if (strlength < 1)
      {
        CCTK_WARN(0, "Need the independent variable names");
      }
      local_indep_names = 
        (const char*)malloc((strlength+1) * sizeof(char));
      ierr = Util_TableGetString(table_handle,
                                 (strlength+1), local_indep_names,
                                 "Independent variable names");
      if (ierr != strlength)
      {
        CCTK_WARN(0, "What happened here?");
      }

      tmpbuffer = EOS_StrSep(&local_indep_names, " ");
      if (tmpbuffer == NULL) /* line is empty */      
      {
        CCTK_WARN(0, "Error in the independent name string; empty???");
      }
      for (var = 0; var < n_entries; ++var)
      {
        while (*tmpbuffer == '\0')
        {
          tmpbuffer = EOS_StrSep(&local_indep_names, " ");
        }
        
        (CallList.Calls[ncall]).indep_vars[var] = NULL;
        for (indep_var = 0; indep_var < n_entries; ++indep_var)
        {
          if ( CCTK_Equals((FnList.Fns[call_number]).indep_names[indep_var],
                           tmpbuffer) )
          {
            (CallList.Calls[ncall]).indep_vars[var] = 
              tmp_indep_arrays[indep_var];
          }
        }
        if ( !((CallList.Calls[ncall]).indep_vars[var]) )
        {
          CCTK_WARN(0, "Failed to set up the independent variables?");
        }
        tmpbuffer = EOS_StrSep(&local_indep_names, " ");
      }

    }
    else
    {
      CCTK_WARN(0, "If setting up a direct array call the independent "
                "arrays\n(key 'Independent arrays') must be set");
    }
    
    
    if (Util_TableQueryValueInfo(table_handle, &type_code, &n_entries,
                                 "Dependent arrays"))
    {
      if (type_code != CCTK_VARIABLE_POINTER)
      {
        CCTK_WARN(0, "The dependent arrays in a direct "
                  "array call must be passed as as a CCTK_POINTER array.");
      }
      if (n_entries > (FnList.Fns[call_number]).n_deps)
      {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "The EOS function '%s' expects at most %d dependent "
                   "variables.\n"
                   "Your direct call list contains %d entries.",
                   fn_name, (FnList.Fns[call_number]).n_deps,
                   n_entries);
      }
      tmp_dep_arrays = (CCTK_REAL**)malloc(n_entries * sizeof(CCTK_REAL*));
      (CallList.Calls[ncall]).dep_vars = 
        (CCTK_REAL**)malloc((FnList.Fns[call_number]).n_deps * 
                            sizeof(CCTK_REAL*));
      (CallList.Calls[ncall]).which_deps_to_set = 
        (CCTK_INT*)malloc((FnList.Fns[call_number]).n_deps * 
                          sizeof(CCTK_INT));
      n_elements = 
        Util_TableGetPointerArray(table_handle,
                                  n_entries,
                                  (CCTK_POINTER*)tmp_dep_arrays,
                                  "Dependent arrays");
      if (n_elements != n_entries)
      {
        CCTK_WARN(0, "What happened here?");
      }
      /* Now we have to order the dependent gfs correctly */
    
      strlength = Util_TableGetString(table_handle,
                                      0, NULL,
                                      "Dependent variable names");
      if (strlength < 1)
      {
        CCTK_WARN(0, "Need the dependent variable names");
      }
      local_dep_names = 
        (const char*)malloc((strlength+1) * sizeof(char));
      ierr = Util_TableGetString(table_handle,
                                 (strlength+1), local_dep_names,
                                 "Dependent variable names");
      if (ierr != strlength)
      {
        CCTK_WARN(0, "What happened here?");
      }

      tmpbuffer = EOS_StrSep(&local_dep_names, " ");
      if (tmpbuffer == NULL) /* line is empty */      
      {
        CCTK_WARN(0, "Error in the independent name string; empty???");
      }
      for (dep_var = 0; dep_var < n_entries; ++var)
      {
        while (*tmpbuffer == '\0')
        {
          tmpbuffer = EOS_StrSep(&local_dep_names, " ");
        }
        
        (CallList.Calls[ncall]).dep_vars[var] = NULL;
        for (dep_var = 0; 
             dep_var < (FnList.Fns[call_number]).n_deps; 
             ++dep_var)
        {
          (CallList.Calls[ncall]).dep_vars[dep_var] = NULL;
          (CallList.Calls[ncall]).which_deps_to_set[dep_var] = 0;
        }
        for (var = 0; var < n_entries; ++var)
        {
          while (*tmpbuffer == '\0')
          {
            tmpbuffer = EOS_StrSep(&local_dep_names, " ");
          }

          for (dep_var = 0; 
               dep_var < (FnList.Fns[call_number]).n_deps; 
               ++dep_var)
          {
            if ( CCTK_Equals((FnList.Fns[call_number]).dep_names[dep_var],
                             tmpbuffer) )
            {
              (CallList.Calls[ncall]).dep_vars[dep_var] = tmp_dep_arrays[var];
              (CallList.Calls[ncall]).which_deps_to_set[dep_var] = 1;
            }
          }
          tmpbuffer = EOS_StrSep(&local_dep_names, " ");
        }
      }
      
    }
    else
    {
      CCTK_WARN(0, "If setting up a direct array call the dependent "
                "arrays\n(key 'Dependent arrays') must be set");
    }
 
  }
  
  
#ifdef EOSBASE_DBG
  printf("An eos call has been setup. The call handle was %d.\n"
         "The eos name is '%s', corresponding to registered eos %d.\n"
         "The registered structure really believes a fn number %d.\n",
         ncall,
         (FnList.Fns)[call_number].eos_name, call_number,
         (CallList.Calls)[ncall].fn_number);
  if ( (CallList.Calls)[ncall].indep_gfs)
  {
    for (var = 0; var < (FnList.Fns)[call_number].n_indeps; ++var)
    {
      printf("Independent gf number %d has index '%d' (name '%s').\n",
             var, (CallList.Calls)[ncall].indep_gfs[var],
             (FnList.Fns)[call_number].indep_names[var]);
    }
    for (var = 0; var < (FnList.Fns)[call_number].n_deps; ++var)
    {
      printf("Dependent gf number %d has index '%d' (to_set = %d, name '%s').\n",
             var, (CallList.Calls)[ncall].dep_gfs[var],
             (CallList.Calls)[ncall].which_deps_to_set[var],
             (FnList.Fns)[call_number].dep_names[var]);
    }
  }
#endif

  return ncall;
}

 /*@@
   @routine    SetGFs
   @date       Fri Jun 10 05:08:20 2005
   @author     Ian Hawke
   @desc 
   The aliased function call by the user thorn that actually sets
   the GFs (in this case).
   The call is passed on to the EOS thorn once the appropriate
   information has been pulled out of the table.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

CCTK_INT SetGFs(const CCTK_POINTER_TO_CONST _GH, const CCTK_INT call_number)
{

  const cGH* const GH = (const cGH* const) _GH;
  
  CCTK_INT ierr = 0;

  CCTK_INT totalsize = 1;
  CCTK_INT dim;
  CCTK_INT fn_number, var;

  for (dim = 0; dim < GH->cctk_dim; ++dim)
  {
    totalsize *= GH->cctk_lsh[dim];
  }
  
  fn_number = (CallList.Calls)[call_number].fn_number;

  for (var = 0; var < (FnList.Fns)[fn_number].n_indeps; ++var)
  {
    (CallList.Calls)[call_number].indep_vars[var] =
      (CCTK_REAL*)CCTK_VarDataPtrI(GH, 0,
                                   (CallList.Calls)[call_number].indep_gfs[var]);
  }
  
  for (var = 0; var < (FnList.Fns)[fn_number].n_deps; ++var)
  {
    if ((CallList.Calls)[call_number].which_deps_to_set[var])
    {
      (CallList.Calls)[call_number].dep_vars[var] =
        (CCTK_REAL*)CCTK_VarDataPtrI(GH, 0,
                                     (CallList.Calls)[call_number].dep_gfs[var]);
    }
    else
    {
      (CallList.Calls)[call_number].dep_vars[var] = NULL;
    }
  }

  ierr = ((FnList.Fns)[fn_number].EOS_fn) 
    ((FnList.Fns)[fn_number].param_table,
     totalsize,
     (CallList.Calls)[call_number].indep_vars,
     (CallList.Calls)[call_number].which_deps_to_set,
     (CallList.Calls)[call_number].dep_vars);

  return ierr;
}

 /*@@
   @routine    SetArrays
   @date       Fri Jun 10 05:09:40 2005
   @author     Ian Hawke
   @desc 
   The aliased function call by the user thorn that actually sets
   the arrays (in this case).
   The call is passed on to the EOS thorn once the appropriate
   information has been pulled out of the table.   
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

CCTK_INT SetArrays(const CCTK_INT call_number)
{
  
  CCTK_INT ierr = 0;

  CCTK_INT fn_number, totalsize;
  
  fn_number = (CallList.Calls)[call_number].fn_number;

  totalsize = (CallList.Calls)[call_number].n_elems;

  ierr = ((FnList.Fns)[fn_number].EOS_fn) 
    ((FnList.Fns)[fn_number].param_table,
     totalsize,
     (CallList.Calls)[call_number].indep_vars,
     (CallList.Calls)[call_number].which_deps_to_set,
     (CallList.Calls)[call_number].dep_vars);
  
  return ierr;
}

CCTK_INT SetInverseArrays(const CCTK_INT call_number, 
                          const CCTK_INT table_handle)
{
  
  CCTK_INT ierr = 0, in_error = 1;
  CCTK_INT fn_number, totalsize, i;
  CCTK_INT unknown_indep_number, known_dep_number;
  CCTK_INT strlength, strlength_recv;
  
  CCTK_REAL guess_factor = 1.01, threshold = 1e-12;
  CCTK_REAL *indep_min, *indep_max, *unknown_indep, 
    *known_result, *tmp_result, *max_result;
  
  char *tmpbuffer;

  fn_number = (CallList.Calls)[call_number].fn_number;

  totalsize = (CallList.Calls)[call_number].n_elems;
  
  /* Find out which independent variable needs setting */

  strlength = Util_TableGetString(table_handle,
                                  0, NULL,
                                  "Unknown independent name");
  if (strlength < 1)
  {
    CCTK_WARN(0, "We need the independent variable name that you"
              "actually want to set\n (key: 'Unknown independent name')");
  }
  tmpbuffer = 
    (const char*)malloc((strlength+1) * sizeof(char));
  strlength_recv = Util_TableGetString(table_handle,
                                       (strlength+1), tmpbuffer,
                                       "Unknown independent name");
  if (strlength_recv != strlength)
  {
    CCTK_WARN(0, "What happened here?");
  }

  unknown_indep_number = -1;
  for (i = 0; i < (FnList.Fns)[fn_number].n_indeps; ++i)
  {
    if (CCTK_Equals(tmpbuffer, (FnList.Fns)[fn_number].indep_names[i]))
    {
      unknown_indep_number = i;
    }
  }
  if (unknown_indep_number < 0)
  {
    CCTK_WARN(0, "Unknown independent variable not recognized");
  }

#ifdef EOSBASE_DBG
  printf("Inverse array call: unknown indep '%s' index %d.\n",
         tmpbuffer, unknown_indep_number);
#endif

  /* Find out which dependent variable is actually known */

  strlength = Util_TableGetString(table_handle,
                                  0, NULL,
                                  "Known dependent name");
  if (strlength < 1)
  {
    CCTK_WARN(0, "We need the dependent variable name that you"
              "have already set\n (key: 'Known dependent name')");
  }
  tmpbuffer = 
    (const char*)malloc((strlength+1) * sizeof(char));
  strlength_recv = Util_TableGetString(table_handle,
                                       (strlength+1), tmpbuffer,
                                       "Known dependent name");
  if (strlength_recv != strlength)
  {
    CCTK_WARN(0, "What happened here?");
  }

  known_dep_number = -1;
  for (i = 0; i < (FnList.Fns)[fn_number].n_deps; ++i)
  {
    if (CCTK_Equals(tmpbuffer, (FnList.Fns)[fn_number].dep_names[i]))
    {
      known_dep_number = i;
    }
  }
  if (known_dep_number < 0)
  {
    CCTK_WARN(0, "Known dependent variable not recognized");
  }

#ifdef EOSBASE_DBG
  printf("Inverse array call: known dep '%s' index %d.\n",
         tmpbuffer, known_dep_number);
#endif

  /* Ensure that the 'known' variable is set */

  (CallList.Calls)[call_number].which_deps_to_set[known_dep_number] = 1;

  /* Juggle... */

  unknown_indep = 
    (CallList.Calls)[call_number].indep_vars[unknown_indep_number];
  known_result = 
    (CallList.Calls)[call_number].dep_vars[known_dep_number];

  indep_min  = (CCTK_REAL*)malloc(totalsize * sizeof(CCTK_REAL));
  indep_max  = (CCTK_REAL*)malloc(totalsize * sizeof(CCTK_REAL));
  tmp_result = (CCTK_REAL*)malloc(totalsize * sizeof(CCTK_REAL));
  max_result = (CCTK_REAL*)malloc(totalsize * sizeof(CCTK_REAL));

  for (i = 0; i < totalsize; ++i)
  {
    indep_min[i] = unknown_indep[i] / guess_factor;
    indep_max[i] = unknown_indep[i] * guess_factor;

#ifdef EOSBASE_DBG
  printf("Inverse array call: index %d (/%d), range [%g, %g].\n",
         i, totalsize, indep_min[i], indep_max[i]);
#endif

  }

  for (i = 0; i < totalsize; ++i)
  {
    unknown_indep[i] = indep_max[i];
  }

  (CallList.Calls)[call_number].dep_vars[known_dep_number] = max_result;
  ierr = ((FnList.Fns)[fn_number].EOS_fn) 
    ((FnList.Fns)[fn_number].param_table,
     totalsize,
     (CallList.Calls)[call_number].indep_vars,
     (CallList.Calls)[call_number].which_deps_to_set,
     (CallList.Calls)[call_number].dep_vars);

  for (i = 0; i < totalsize; ++i)
  {
    unknown_indep[i] = indep_min[i];
  }

  (CallList.Calls)[call_number].dep_vars[known_dep_number] = tmp_result;
  ierr = ((FnList.Fns)[fn_number].EOS_fn) 
    ((FnList.Fns)[fn_number].param_table,
     totalsize,
     (CallList.Calls)[call_number].indep_vars,
     (CallList.Calls)[call_number].which_deps_to_set,
     (CallList.Calls)[call_number].dep_vars);

  for (i = 0; i < totalsize; ++i)
  {
    if ( ((known_result[i] - tmp_result[i]) * 
          (known_result[i] - max_result[i])) > 0 )
    {
#ifdef EOSBASE_DBG
  printf("Inverse array call: index %d (/%d), bracketing failed?\n"
         "Range [%g, %g], result %g.\n",
         i, totalsize, tmp_result[i], max_result[i], known_result[i]);
#endif

      CCTK_WARN(0, "Failed to bracket root!");
    }
    unknown_indep[i] = 0.5 * (indep_min[i] + indep_max[i]);
  }
  

  while (in_error)
  {

    ierr = ((FnList.Fns)[fn_number].EOS_fn) 
      ((FnList.Fns)[fn_number].param_table,
       totalsize,
       (CallList.Calls)[call_number].indep_vars,
       (CallList.Calls)[call_number].which_deps_to_set,
       (CallList.Calls)[call_number].dep_vars);

    in_error = 0;

    for (i = 0; i < totalsize; ++i)
    {
      if ( ((known_result[i] - tmp_result[i]) * 
            (known_result[i] - max_result[i])) < 0 )
      {
        indep_min[i] = unknown_indep[i];
      }
      else
      {
        indep_max[i] = unknown_indep[i];
        max_result[i] = tmp_result[i];
      }
      unknown_indep[i] = 0.5 * (indep_min[i] + indep_max[i]);
      if (indep_max[i] - indep_min[i] > threshold)
      {
        in_error = 1;
      }
    }
  
  }
  
  (CallList.Calls)[call_number].dep_vars[known_dep_number] = known_result;

  free(indep_min); indep_min = NULL;
  free(indep_max); indep_max = NULL;
  free(tmp_result); tmp_result = NULL;
  free(max_result); max_result = NULL;

  return ierr;
}
