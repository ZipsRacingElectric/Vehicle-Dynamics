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
** Simple Brake Model for a four wheel vehicle
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    Brake_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

#define NWHEEL 4

static const char ThisModelClass[] = "Brake";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


struct tMyModel {
    double	TrqDistrib[NWHEEL];	/* Model Parameter */
    double	Trq_WB[NWHEEL];		/* Brake torque [Nm] */
};



/*
** MyModel_New ()
**
** Initialising the model
**
** Call:
** - one times at the beginning of every TestRun
*/

static void *
MyModel_New (tInfos *Inf, tBrakeCfgIF *CfgIF, const char *KindKey)
{
    struct tMyModel *mp = NULL;
    const char  *ModelKind, *key;
    char	MsgPre[64];
    double	dvec[NWHEEL];
    int		i, VersionId = 0;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);
    if (CfgIF->nWheels != NWHEEL) {
    	LogErrF (EC_Init, "%s: wheel number mismatch (brake %d, vehicle %d)", MsgPre,
		 NWHEEL,CfgIF->nWheels);
    	return NULL;
    }

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_Brake, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*)calloc(1, sizeof(struct tMyModel));

    /*** Get parameters of MyModel  */
    key = "Brake.Pedal2Trq";
    iGetTable(Inf, key, dvec, NWHEEL, 1, &i);
    if (i != NWHEEL) {
	LogErrF (EC_Init, "%s: Unsupported argument for '%s'", MsgPre, key);
	goto ErrorReturn;
    }

    for (i=0; i<NWHEEL; i++)
	mp->TrqDistrib[i] = dvec[i];

    return mp;

    ErrorReturn:
	free (mp);
	return NULL;
}


static void
MyModel_DeclQuants_dyn (struct tMyModel *mp, int park)
{
    int i;
    char buf[64];
    static struct tMyModel MyModel_Dummy = {{0}};
    if (park)
	mp = &MyModel_Dummy;

    /* Define here dict entries for dynamically allocated variables. */

    for (i=0; i<NWHEEL; i++) {
	sprintf (buf, "MyModel.Trq_%s", Vehicle_TireNo_Str(i));
	DDefDouble4 (NULL, buf, "Nm", &mp->Trq_WB[i], DVA_None);
    }
}


/*
** MyModel_DeclQuants ()
**
** Defining DataDictionary Quantities
**
** Call:
** - at the beginning of every TestRun, after the call for the _New function
*/

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


/*
** MyModel_Calc ()
**
** Calculation of the model
**
** Call:
** - every cycle
*/


static int
MyModel_Calc (
    void	    *MP,	/* model parameter handle	*/
    struct tBrakeIF *IF,	/* brake interface structure	*/
    double	    dt)		/* time step			*/

{
    struct tMyModel *mp = (struct tMyModel *)MP;
    int i;

    for (i=0; i < NWHEEL; i++) {
        mp->Trq_WB[i] = mp->TrqDistrib[i] * IF->Pedal;
        IF->Trq_WB[i] = mp->Trq_WB[i];
    }

    return 0;
}


/*
** MyModel_Delete ()
**
** Uninitialising the model
**
** Call:
** - at the end of every TestRun
*/

static void
MyModel_Delete (void *MP)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Park the dict entries for dynamically allocated variables before deleting */
    MyModel_DeclQuants_dyn (mp, 1);
    free (mp);
}



/*
** Brake_Register_MyModel
**
*/

int
Brake_Register_MyModel (void)
{
    tModelClassDescr	mp;

    memset (&mp, 0, sizeof(mp));
    mp.Brake.New = 	   MyModel_New;
    mp.Brake.Calc = 	   MyModel_Calc;
    mp.Brake.Delete = 	   MyModel_Delete;
    mp.Brake.DeclQuants =  MyModel_DeclQuants;
    mp.Brake.VersionId =   ThisVersionId;
    /* Should only be used if the model doesn't read params from extra files */
    mp.Brake.ParamsChanged = ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_Brake, ThisModelKind, &mp);
}
