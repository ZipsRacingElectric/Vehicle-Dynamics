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
** Simple external suspension Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    int SuspExtFrcs_Register_MyModel (void);
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Susp.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "SuspExtFrcs";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;

#define N_SUSPENSIONS 4

/* Modellparameter (statische) */
struct tMyModel {
    struct tMyFrcSpring {
	double len0;
	double Frc0;
	double dFrc_dl;
    } Spring[N_SUSPENSIONS];

    struct tMyFrcDamp {
	double dFrc_dlp;
    } Damp_Push[N_SUSPENSIONS], Damp_Pull[N_SUSPENSIONS];

    struct tMyFrcBuffer {
    	double len0;
	double dFrc_dl;
    } Buf_Push[N_SUSPENSIONS], Buf_Pull[N_SUSPENSIONS];

    struct tMyFrcStabi {
	double dFrc_dl;
    } Stabi[N_SUSPENSIONS];
};


static void
MyModel_DeclQuants_dyn (struct tMyModel *mp, int park)
{
    static struct tMyModel MyModel_Dummy;
    memset (&MyModel_Dummy, 0, sizeof(struct tMyModel));
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

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_SuspExtFrcs, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    if (CfgIF->nWheels != N_SUSPENSIONS) {
	LogErrF (EC_Init, "%s: supports only a vehicle with 4 wheels", MsgPre);
	return NULL;
    }

    mp = (struct tMyModel*)calloc(1, sizeof(*mp));

    for (iS=0; iS < N_SUSPENSIONS; iS++) {
	char Key[32];
	const char *s = Vehicle_TireNo_Str(iS);
	double	*dv = NULL;

	/* key = <len0> <Frc0> <dFrc/dl> */
	sprintf(Key, "SFH.Spring%s", s);
	dv = iGetTable2(Inf, Key, 3, &n);
	if (dv != NULL) {
	    mp->Spring[iS].len0 =	dv[0];
	    mp->Spring[iS].Frc0 =	dv[1];
	    mp->Spring[iS].dFrc_dl =	dv[2];
	    free(dv);
	} else {
	    LogErrF (EC_Init, "%s: Unsupported argument for '%s'", MsgPre, Key);
	    goto ErrorReturn;
	}

	/* key = <dF/dlp> */
	sprintf(Key, "SFH.Damp_Push%s", s);
	mp->Damp_Push[iS].dFrc_dlp = iGetDbl(Inf, Key);

	sprintf(Key, "SFH.Damp_Pull%s", s);
	mp->Damp_Pull[iS].dFrc_dlp = iGetDbl(Inf, Key);
	
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

	/* key = <dF/dl> */
	sprintf(Key, "SFH.Stabi%s", s);
	mp->Stabi[iS].dFrc_dl = iGetDbl(Inf, Key);
    }

    return mp;

    ErrorReturn:
	free(mp);
    	return NULL;
}


static int
MyModel_Calc (void *MP, tSuspExtFrcsIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    int	iS;

    for (iS=0; iS < N_SUSPENSIONS; iS++) {
	IF->FrcSpring[iS]
	    = mp->Spring[iS].dFrc_dl * (mp->Spring[iS].len0 - IF->lSpring[iS])
	    + mp->Spring[iS].Frc0;

	if (IF->vDamp[iS] >= 0.0) {
	    IF->FrcDamp[iS] = -mp->Damp_Pull[iS].dFrc_dlp * IF->vDamp[iS];
	} else {
	    IF->FrcDamp[iS] = -mp->Damp_Push[iS].dFrc_dlp * IF->vDamp[iS];
	}

	IF->FrcBuf[iS] = 0.0;
	if (IF->lBuf[iS] < mp->Buf_Push[iS].len0) {
	    IF->FrcBuf[iS] = mp->Buf_Push[iS].dFrc_dl * (mp->Buf_Push[iS].len0 - IF->lBuf[iS]);
	} else if (IF->lBuf[iS] > mp->Buf_Pull[iS].len0) {
	    IF->FrcBuf[iS] = - mp->Buf_Pull[iS].dFrc_dl * (IF->lBuf[iS] - mp->Buf_Pull[iS].len0);
	}
    }

    for (iS=0; iS < N_SUSPENSIONS; iS += 2) {
	/* Frc = c * (right - left) */
	IF->FrcStabi[iS]   = mp->Stabi[iS].dFrc_dl * (IF->lStabi[iS+1] - IF->lStabi[iS]);
	IF->FrcStabi[iS+1] = -IF->FrcStabi[iS];
    }

    return 0;
}


static void
MyModel_Delete (void *MP)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Park the dict entries for dynamically allocated variables before deleting */
    MyModel_DeclQuants_dyn (mp, 1);
    free (mp);
}


int
SuspExtFrcs_Register_MyModel (void)
{
    tModelClassDescr m;

    memset (&m, 0, sizeof(m));
    m.SuspExtFrcs.VersionId =		ThisVersionId;
    m.SuspExtFrcs.New =			MyModel_New;
    m.SuspExtFrcs.Calc =		MyModel_Calc;
    m.SuspExtFrcs.Delete =		MyModel_Delete;
    m.SuspExtFrcs.DeclQuants =		MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.SuspExtFrcs.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_SuspExtFrcs, ThisModelKind, &m);
}
