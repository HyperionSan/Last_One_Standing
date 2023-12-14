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

 /*@@
   @routine    EOS_StrSep
   @date       Wed Mar  2 14:57:50 2005
   @author     Ian Hawke
   @desc 
   A wrapper to Util_StrSep that returns the last token if there 
   are no more delimiters. Note that the documentation is incorrect;
   delim_set (here delim) cannot be a set but must be a single 
   character.
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/


const char * EOS_StrSep(const char ** string_ptr, const char * delim);

const char * EOS_StrSep(const char ** string_ptr, const char * delim)
{
  
  const char * token;
  const char ** string;
  
  string = string_ptr;
  
  if (string == NULL)
  {
    return NULL;
  }
  
  token = Util_StrSep(string_ptr, delim);
  
  if ( (token == NULL) && (strlen(*string) > 0) && 
       (!(CCTK_Equals(*string, "\n"))) )
  {
    string_ptr = string;
    token = Util_StrSep(string_ptr, "\n");
    if ( (token == NULL) && (strlen(*string) > 0) )
    {
      string_ptr = '\0';
      token = *string;
    } 
  }

#ifdef EOS_GTR_DEBUG
  printf("EOS_StrSep: token '%s' and string_ptr '%s' (%s)\n",
         token, *string_ptr, *string);
#endif

  return token;
}
