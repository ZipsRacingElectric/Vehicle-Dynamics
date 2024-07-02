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
** Simple suspension Model for damper 
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    int Susp_Damper_Register_MyModel (void);
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

static const char ThisModelClass[] = "Susp_Damper";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;

/* model parameters */
struct tMyModel {
    struct tMyFrcDamp {
	double dFrc_dlp;
    } Damp_Push, Damp_Pull;
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
MyModel_DeclQuants (void *MP, int SuspModID)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    if (mp == NULL) {
	/* Define here dict entries for non-dynamically allocated (static) variables. */

    } else {
	MyModel_DeclQuants_dyn (mp, 0);
    }
}


static void *
MyModel_New (struct tInfos *Inf, struct tSusp_DamperCfgIF *CfgIF, const char *KindKey)
{
    struct tMyModel *mp = NULL;
    int VersionId = 0;
    char MsgPre[64];
    const char *ModelKind;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_Susp_Damper, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    mp = (struct tMyModel*)calloc(1, sizeof(*mp));

    char Key[32];
    const char *s = SuspMod_SuspPre(CfgIF->SuspModID);

    /* key = <dF/dlp> */
    sprintf(Key, "SFH.%s.Damp_Push", s);
    mp->Damp_Push.dFrc_dlp = iGetDbl(Inf, Key);

    sprintf(Key, "SFH.%s.Damp_Pull", s);
    mp->Damp_Pull.dFrc_dlp = iGetDbl(Inf, Key);

    Susp_Damper_GetCfgOutIF(Inf, CfgIF, ModelKind);

    return mp;
}


static int
MyModel_Calc (void *MP, tSusp_DamperIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel*)MP;

    if (IF->vel >= 0.0) {
	IF->Force = -mp->Damp_Pull.dFrc_dlp * IF->vel;
    } else {
	IF->Force = -mp->Damp_Push.dFrc_dlp * IF->vel;
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
Susp_Damper_Register_MyModel (void)
{
    tModelClassDescr m;

    memset (&m, 0, sizeof(m));
    m.Susp_Damper.VersionId =	  ThisVersionId;
    m.Susp_Damper.New =		  MyModel_New;
    m.Susp_Damper.Calc =	  MyModel_Calc;
    m.Susp_Damper.Delete =	  MyModel_Delete;
    m.Susp_Damper.DeclQuants =	  MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.Susp_Damper.ParamsChanged = ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_Susp_Damper, ThisModelKind, &m);
}
