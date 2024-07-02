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
** Simple VehicleControl Model to demonstrate the manipulation of Driver Gas
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    VehicleControl_Register_MyModel ()
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

struct tMyModel {

    /* parameters of MyVehicleControl */
    double GasFac;

    /* variables of MyVehicleControl */
    double Gas;
};

static const char ThisModelClass[] = "VehicleControl";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


/******************************************************************************/


static void
MyModel_DeclQuants_dyn (struct tMyModel *mp, int park)
{
    static struct tMyModel MyModel_Dummy = {0};
    if (park)
	mp = &MyModel_Dummy;

    /* Define here dict entries for dynamically allocated variables. */

    DDefDouble4 (NULL, "MyVehicleControl.Gas", "", &mp->Gas, DVA_None);
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


static int
MyModel_Calc (void *MP, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Manipulation of VehicleControl */
    mp->Gas = mp->GasFac * VehicleControl.Gas;
    VehicleControl.Gas = mp->Gas;

    return 0;
}


static void *
MyModel_New (struct tInfos *Inf, const char *KindKey)
{
    struct tMyModel *mp = NULL;
    const char *ModelKind, *key;
    char MsgPre[64];
    int VersionId = 0;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_VehicleControl, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel *)calloc(1, sizeof(*mp));

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    key = "MyVC.GasFac";
    mp->GasFac = iGetDbl(Inf, key);
    if (mp->GasFac < 0.0 || mp->GasFac > 1.0) {
	LogErrF(EC_Init, "%s: gas factor '%s' must be 0...1", MsgPre, key);
	goto ErrorReturn;
    }

    return mp;

    ErrorReturn:
	free (mp);
	return NULL;
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
VehicleControl_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.VehicleControl.VersionId =	ThisVersionId;	
    m.VehicleControl.New    = 		MyModel_New;
    m.VehicleControl.Calc   = 		MyModel_Calc;
    m.VehicleControl.DeclQuants   = 	MyModel_DeclQuants;
    m.VehicleControl.Delete = 		MyModel_Delete;
    /* Should only be used if the model doesn't read params from extra files */
    m.VehicleControl.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_VehicleControl, ThisModelKind, &m);
}
