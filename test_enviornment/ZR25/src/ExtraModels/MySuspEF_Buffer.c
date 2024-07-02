/*
******************************************************************************
**  CarMaker - Version 12.0.1
**  Vehicle Dynamics Simulation Toolkit
**
**  Copyright (C)   IPG Automotive GmbH
**                  Bannwaldallee 60             Phone  +49.721.98520.0
**                  76185 Karlsruhe              Fax    +49.721.98520.99
**                  Germany                      WWW    www.ipg-automotive.com
******************************************************************************
**
** Simple suspension Model for external buffer
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    int SuspEF_Buffer_Register_MyModel (void);
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Car.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"
#include "Log.h"

static const char ThisModelClass[] = "SuspEF_Buffer";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;

#define N_SUSPENSIONS 4

/* Modellparameter (statische) */
struct tMyModel {
    struct tMyFrcBuffer {
    	double len0;
	double dFrc_dl;
    } Buf_Push[N_SUSPENSIONS], Buf_Pull[N_SUSPENSIONS];

};


static void
MyModel_DeclQuants_dyn (struct tMyModel *mp, int park)
{
    static struct tMyModel MyModel_Dummy = { { { 0, 0 } } };
    if (park)
	mp = &MyModel_Dummy;

    /* Define here dict entries for dynamically allocated variables. */
}


static void
MyModel_DeclQuants (void *MP)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    if (mp == NULL) {
	/* Define here dict entries for non-dynamically allocated (static) variables. */

    } else {
	MyModel_DeclQuants_dyn (mp, 0);
    }
}


static void *
MyModel_New (struct tInfos *Inf, struct tSuspExtFrcsCfgIF *CfgIF, const char *KindKey)
{
    struct tMyModel *mp = NULL;
    int iS, n, VersionId = 0;
    char MsgPre[64];
    const char *ModelKind;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_SuspEF_Buffer, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    mp = (struct tMyModel*)calloc(1, sizeof(*mp));

    for (iS=0; iS < N_SUSPENSIONS; iS++) {
	char Key[32];
	const char *s = Vehicle_TireNo_Str(iS);
	double	*dv = NULL;

	/* key = <len0>  <dFrc/dl> */
	sprintf(Key, "SFH.Buf_Push%s", s);
	dv = iGetTable2(Inf, Key, 2, &n);
	if (dv != NULL) {
	    mp->Buf_Push[iS].len0 =	dv[0];
	    mp->Buf_Push[iS].dFrc_dl =	dv[1];
	    free(dv);
	} else {
	    LogErrF (EC_Init, "%s: Unsupported argument for '%s'", MsgPre, Key);
	    goto ErrorReturn;
	}
	
	sprintf(Key, "SFH.Buf_Pull%s", s);
	dv = iGetTable2(Inf, Key, 2, &n);
	if (dv != NULL) {
	    mp->Buf_Pull[iS].len0 =	dv[0];
	    mp->Buf_Pull[iS].dFrc_dl =	dv[1];
	    free(dv);
	} else {
	    LogErrF (EC_Init, "%s: Unsupported argument for '%s'", MsgPre, Key);
	    goto ErrorReturn;
	}
    }

    return mp;

    ErrorReturn:
	if (mp != NULL)
	    free(mp);
    	return NULL;
}


static int
MyModel_Calc (void *MP, tSuspEF_BufferIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel*)MP;
    int	iS;

    for (iS=0; iS < N_SUSPENSIONS; iS++) {
	IF->FrcBuf[iS] = 0.0;
	if (IF->lBuf[iS] < mp->Buf_Push[iS].len0) {
	    IF->FrcBuf[iS] = mp->Buf_Push[iS].dFrc_dl * (mp->Buf_Push[iS].len0 - IF->lBuf[iS]);
	} else if (IF->lBuf[iS] > mp->Buf_Pull[iS].len0) {
	    IF->FrcBuf[iS] = - mp->Buf_Pull[iS].dFrc_dl * (IF->lBuf[iS] - mp->Buf_Pull[iS].len0);
	}
    }

    return 0;
}


static void
MyModel_Delete (void *MP)
{
    struct tMyModel *mp = (struct tMyModel*)MP;

    if (mp != NULL)
	free (mp);
    mp = NULL;
}


int
SuspEF_Buffer_Register_MyModel (void)
{
    tModelClassDescr m;

    memset (&m, 0, sizeof(m));
    m.SuspEF_Buffer.VersionId =		ThisVersionId;
    m.SuspEF_Buffer.New =		MyModel_New;
    m.SuspEF_Buffer.Calc =		MyModel_Calc;
    m.SuspEF_Buffer.Delete =		MyModel_Delete;
    m.SuspEF_Buffer.DeclQuants =	MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.SuspEF_Buffer.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_SuspEF_Buffer, ThisModelKind, &m);
}
