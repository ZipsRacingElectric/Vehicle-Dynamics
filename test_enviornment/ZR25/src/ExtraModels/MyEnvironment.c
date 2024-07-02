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
** Simple environment Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register())
**
**    Environment_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelKind[] = "MyEnvironment";
static const int  ThisVersionId   = 1;


struct tMyModel {
    /* Parameters */
    double	Temp;		/* Air temperature [K]   */
    double	Pressure;	/* Air pressure    [bar] */
};


static void
MyModel_DeclQuants_dyn (struct tMyModel *mp, int park)
{
    static struct tMyModel MyModel_Dummy = {0};
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
MyModel_New (tInfos *Inf, const char *KindKey)
{
    struct tMyModel *mp = NULL;
    const char *ModelKind;
    int VersionId = 0;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_Environment, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*)calloc(1, sizeof(*mp));

    mp->Temp     = iGetDblOpt (Inf, "Env.MyEnvironment.Temp",     293.15);
    mp->Pressure = iGetDblOpt (Inf, "Env.MyEnvironment.Pressure", 1.013);

    return mp;
}


static int
MyModel_Calc (void *MP, tEnvironmentIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    IF->Temperature = mp->Temp;
    IF->AirPressure = mp->Pressure;

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
Environment_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.Environment.VersionId =	ThisVersionId;
    m.Environment.New =		MyModel_New;
    m.Environment.Calc =	MyModel_Calc;
    m.Environment.Delete =	MyModel_Delete;
    m.Environment.DeclQuants =	MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.Environment.ParamsChanged = ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_Environment, ThisModelKind, &m);
}
