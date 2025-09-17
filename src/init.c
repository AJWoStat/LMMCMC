#include <stdlib.h>
#include <R_ext/Rdynload.h>

#include "lm_mcmc_function.h"

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef R_CallDef[] = {
  CALLDEF(lm_mcmc_function, 17),
  {NULL, NULL, 0}
};

void R_init_splines(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, R_CallDef, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
