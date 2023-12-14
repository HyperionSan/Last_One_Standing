/* -------------------------------------------------------------------------- */
/* UMFPACK Version 4.0 (Apr 11, 2002), Copyright (c) 2002 by Timothy A.       */
/* Davis.  All Rights Reserved.  See README for License.                      */
/* email: davis@cise.ufl.edu    CISE Department, Univ. of Florida.            */
/* web: http://www.cise.ufl.edu/research/sparse/umfpack                       */
/* -------------------------------------------------------------------------- */

#if defined (UMF_MALLOC_COUNT) || !defined (NDEBUG)
GLOBAL extern Int UMF_malloc_count ;
#endif

GLOBAL void *UMF_malloc
(
    Int n_objects,
    size_t size_of_object
) ;

