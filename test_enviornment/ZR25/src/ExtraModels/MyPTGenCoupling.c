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
** Simple powertrain Generic coupling Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    PTGenCoupling_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "PT Coupling";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


struct tMyModel {
    /* Parameters */
    tPTGen_CplDiffPos	CplPos;
    double		k;	/* Friction coefficient */
};


static void
MyModel_DeclQuants_dyn (struct tMyModel *mp, int park)
{
    static struct tMyModel MyModel_Dummy = { PTGen_Pos_Unknown };
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

/* Model output parameters in the configuration struct CfgIF, which are required
   by CarMaker, are read in before the MyModel_New() function.
   - The parametrization of these parameters is supported by the GUI.
   - These output parameters can be used internally by the model in same way like
     the input parameters
*/
static void *
MyModel_New (
    struct tInfos		*Inf,
    struct tPTGenCouplingCfgIF	*CfgIF,
    const char      		*KindKey,
    const char			*Pre)
{
    struct tMyModel *mp = NULL;
    const char *ModelKind, *key;
    char MsgPre[64];
    int VersionId = 0;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_PTGenCoupling, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*)calloc(1,sizeof(*mp));
    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    /* get CfgIF parameters */
    if (Coupling_GetCfgOutIF (Inf, CfgIF, ModelKind, Pre) != 0)
	goto ErrorReturn;

    /* CfgIF -> Model */
    mp->CplPos = CfgIF->CplPos;

    key = "MyPTGenCoupling.rotv2Trq";
    mp->k = iGetDbl(Inf, key);
    if (mp->k<0) {
	LogErrF (EC_Init, "%s: parameter '%s' must be positive or zero",
	    MsgPre, key);
	goto ErrorReturn;
    }

    /* CfgIF output: verification if the parametrization corresponds to the model */
    if (CfgIF->CplType != PTGen_CplType_NotLockable) {
	LogErrF (EC_Init, "%s: model supports only a not lockable coupling type", MsgPre);
	goto ErrorReturn;
    }

    return mp;

    ErrorReturn:
	free (mp);
	return NULL;
}


static int
MyModel_Calc (void *MP, struct tPTGenCouplingIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    IF->Trq_A2B = mp->k * IF->drotv_A2B;

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
PTGenCoupling_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.PTGenCoupling.VersionId =		ThisVersionId;
    m.PTGenCoupling.New =		MyModel_New;
    m.PTGenCoupling.Calc =		MyModel_Calc;
    m.PTGenCoupling.Delete =		MyModel_Delete;
    m.PTGenCoupling.DeclQuants =	MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.PTGenCoupling.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_PTGenCoupling, ThisModelKind, &m);
}
