#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _pseudoCure_gee(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _pseudoCure_geeCV(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _pseudoCure_pgee(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _pseudoCure_pseudoKM(SEXP, SEXP);
extern SEXP _pseudoCure_pseudoKM1(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_pseudoCure_gee",       (DL_FUNC) &_pseudoCure_gee,        8},
    {"_pseudoCure_geeCV",     (DL_FUNC) &_pseudoCure_geeCV,     13},
    {"_pseudoCure_pgee",      (DL_FUNC) &_pseudoCure_pgee,      12},
    {"_pseudoCure_pseudoKM",  (DL_FUNC) &_pseudoCure_pseudoKM,   2},
    {"_pseudoCure_pseudoKM1", (DL_FUNC) &_pseudoCure_pseudoKM1,  3},
    {NULL, NULL, 0}
};

void R_init_pseudoCure(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
