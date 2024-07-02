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
** Simple suspension Model for external damper 
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    int SuspEF_Damper_Register_MyModel (void);
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

static const char ThisModelClass[] = "SuspEF_Damper";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;

#define N_SUSPENSIONS 4

/* Modellparameter (statische) */
struct tMyModel {
    struct tMyFrcDamp {
	double dFrc_dlp;
    } Damp_Push[N_SUSPENSIONS], Damp_Pull[N_SUSPENSIONS];
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
    int iS, VersionId = 0;
    char MsgPre[64];
    const char *ModelKind;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_SuspEF_Damper, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    mp = (struct tMyModel*)calloc(1, sizeof(*mp));

    for (iS=0; iS < N_SUSPENSIONS; iS++) {
	char Key[32];
	const char *s = Vehicle_TireNo_Str(iS);

	/* key = <dF/dlp> */
	sprintf(Key, "SFH.Damp_Push%s", s);
	mp->Damp_Push[iS].dFrc_dlp = iGetDbl(Inf, Key);

	sprintf(Key, "SFH.Damp_Pull%s", s);
	mp->Damp_Pull[iS].dFrc_dlp = iGetDbl(Inf, Key);
    }

    return mp;
}


static int
MyModel_Calc (void *MP, tSuspEF_DamperIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel*)MP;
    int	iS;

    for (iS=0; iS < N_SUSPENSIONS; iS++) {
	if (IF->vDamp[iS] >= 0.0) {
	    IF->FrcDamp[iS] = -mp->Damp_Pull[iS].dFrc_dlp * IF->vDamp[iS];
	} else {
	    IF->FrcDamp[iS] = -mp->Damp_Push[iS].dFrc_dlp * IF->vDamp[iS];
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
SuspEF_Damper_Register_MyModel (void)
{
    tModelClassDescr m;

    memset (&m, 0, sizeof(m));
    m.SuspEF_Damper.VersionId =		ThisVersionId;
    m.SuspEF_Damper.New =		MyModel_New;
    m.SuspEF_Damper.Calc =		MyModel_Calc;
    m.SuspEF_Damper.Delete =		MyModel_Delete;
    m.SuspEF_Damper.DeclQuants =	MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.SuspEF_Damper.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_SuspEF_Damper, ThisModelKind, &m);
}
