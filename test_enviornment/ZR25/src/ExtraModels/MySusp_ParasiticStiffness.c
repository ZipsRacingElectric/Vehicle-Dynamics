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
** Simple suspension Model parasitic stiffness
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    int Susp_ParasiticStiffness_Register_MyModel (void);
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

static const char ThisModelClass[] = "Susp_ParasiticEffects";
static const char ThisModelKind[]  = "MyStiffnessModel";
static const int  ThisVersionId    = 1;


/* model parameters */
struct tMyModel {
    struct tMyFrcParasiticStiffness{
	double dFrc_dl;
    } ParasiticStiffness;
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
MyModel_New (struct tInfos *Inf, struct tSusp_ParasiticEffectsCfgIF *CfgIF, const char *KindKey)
{
    struct tMyModel *mp = NULL;
    int VersionId = 0, rv = 0;
    char MsgPre[64];
    const char *ModelKind;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_Susp_ParasiticEffects, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    mp = (struct tMyModel*)calloc(1, sizeof(*mp));

    char Key[32];
    const char *s = SuspMod_SuspPre(CfgIF->SuspModID);

    /* key = <dF/dl> */
    sprintf(Key, "SFH.%s.Parasitic.Stiffness", s);
    mp->ParasiticStiffness.dFrc_dl = iGetDbl(Inf, Key);

    rv = Susp_ParasiticStiffness_GetCfgOutIF (Inf, CfgIF, ModelKind);
    if (rv < 0)
        goto ErrorReturn;

    return mp;

    ErrorReturn:
    if (mp != NULL)
	free(mp);
    return NULL;
}


static int
MyModel_Calc (void *MP, tSusp_ParasiticEffectsIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel*)MP;

    IF->Force   = mp->ParasiticStiffness.dFrc_dl * IF->tWheelCarrier_z;

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
Susp_ParasiticStiffness_Register_MyModel (void)
{
    tModelClassDescr m;

    memset (&m, 0, sizeof(m));
    m.Susp_ParasiticEffects.VersionId =	    	ThisVersionId;
    m.Susp_ParasiticEffects.New =		MyModel_New;
    m.Susp_ParasiticEffects.Calc =		MyModel_Calc;
    m.Susp_ParasiticEffects.Delete =		MyModel_Delete;
    m.Susp_ParasiticEffects.DeclQuants =	MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.Susp_ParasiticEffects.ParamsChanged = ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_Susp_ParasiticEffects, ThisModelKind, &m);
}
