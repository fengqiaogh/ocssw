#ifndef  _L1_OSMI_H
#define  _L1_OSMI_H

#include <stdint.h>
#include "l1.h"

int closel1_osmi(filehandle *l1file);
int openl1_osmi(filehandle *l1file);
int readl1_osmi(filehandle *l1file, int32_t recnum, l1str *l1rec);

#endif
