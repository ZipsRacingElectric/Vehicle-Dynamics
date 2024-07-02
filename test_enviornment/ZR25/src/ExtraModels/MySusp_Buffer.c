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
** Simple suspension Model for buffer
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    int Susp_Buffer_Register_MyModel (void);
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

static const char ThisModelClass[] = "Susp_Buffer";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


/* model parameters */
struct tMyModel {
    struct tMyFrcBuffer {
    	double len0;
	double dFrc_dl;
    } Buf_Push, Buf_Pull;

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
MyModel_New (struct tInfos *Inf, struct tSusp_BufferCfgIF *CfgIF, const char *KindKey)
{
    struct tMyModel *mp = NULL;
    int n, VersionId = 0;
    char MsgPre[64];
    const char *ModelKind;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_Susp_Buffer, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    mp = (struct tMyModel*)calloc(1, sizeof(*mp));

    char Key[32];
    const char *s = SuspMod_SuspPre(CfgIF->SuspModID);
    double	*dv = NULL;

    /* key = <len0>  <dFrc/dl> */
    sprintf(Key, "SFH.%s.Buf_Push", s);
    dv = iGetTable2(Inf, Key, 2, &n);
    if (dv != NULL) {
	mp->Buf_Push.len0 =	dv[0];
	mp->Buf_Push.dFrc_dl =	dv[1];
	free(dv);
    } else {
	LogErrF (EC_Init, "%s: Unsupported argument for '%s'", MsgPre, Key);
	goto ErrorReturn;
    }
    
    sprintf(Key, "SFH.%s.Buf_Pull", s);
    dv = iGetTable2(Inf, Key, 2, &n);
    if (dv != NULL) {
	mp->Buf_Pull.len0 =	dv[0];
	mp->Buf_Pull.dFrc_dl =	dv[1];
	free(dv);
    } else {
	LogErrF (EC_Init, "%s: Unsupported argument for '%s'", MsgPre, Key);
	goto ErrorReturn;
    }

    Susp_Buffer_GetCfgOutIF (Inf, CfgIF, ModelKind);

    return mp;

    ErrorReturn:
	if (mp != NULL)
	    free(mp);
    	return NULL;
}


static int
MyModel_Calc (void *MP, tSusp_BufferIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel*)MP;

    IF->Force = 0.0;
    if (IF->length < mp->Buf_Push.len0) {
	IF->Force = mp->Buf_Push.dFrc_dl * (mp->Buf_Push.len0 -IF->length);
    } else if (IF->length > mp->Buf_Pull.len0) {
	IF->Force = mp->Buf_Pull.dFrc_dl * (mp->Buf_Pull.len0 -IF->length);
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
Susp_Buffer_Register_MyModel (void)
{
    tModelClassDescr m;

    memset (&m, 0, sizeof(m));
    m.Susp_Buffer.VersionId =	  ThisVersionId;
    m.Susp_Buffer.New =		  MyModel_New;
    m.Susp_Buffer.Calc =	  MyModel_Calc;
    m.Susp_Buffer.Delete =	  MyModel_Delete;
    m.Susp_Buffer.DeclQuants =	  MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.Susp_Buffer.ParamsChanged = ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_Susp_Buffer, ThisModelKind, &m);
}
