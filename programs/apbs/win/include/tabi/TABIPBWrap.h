#ifndef H_TABIPB_WRAP_H
#define H_TABIPB_WRAP_H

#include "TABIPBStruct.h"
#include "generic/valist.h"

#ifdef __cplusplus
extern "C" {
#endif

struct TABIPBOutput runTABIPBWrapAPBS(struct TABIPBInput tabipbIn, Valist* APBSMolecule);

#ifdef __cplusplus
}
#endif

#endif
