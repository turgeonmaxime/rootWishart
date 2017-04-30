#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
 Check these declarations against the C/Fortran source code.
 */

/* .Call calls */
extern SEXP rootWishart_doubleWishart_raw(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP rootWishart_singleWishart_raw(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"rootWishart_doubleWishart_raw", (DL_FUNC) &rootWishart_doubleWishart_raw, 5},
    {"rootWishart_singleWishart_raw", (DL_FUNC) &rootWishart_singleWishart_raw, 4},
    {NULL, NULL, 0}
};

void R_init_rootWishart(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
