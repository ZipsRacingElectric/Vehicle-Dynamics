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
** Simple Hydraulic Brake Controller Model for a four wheel vehicle
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    HydBrakeCU_Register_MyModel ();
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

static const char ThisModelClass[] = "Brake.Control";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;

struct tMyModel {
    double	TrqDistrib[NWHEEL];	/* Model Parameter */
    double	dTrqSum[NWHEEL];	/* Brake torque [Nm] */
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
MyModel_New (tInfos *Inf, tHydBrakeCU_CfgIF *CfgIF, const char *KindKey)
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

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_HydBrakeControl, KindKey,
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
    static struct tMyModel MyModel_Dummy = {{0}};
    if (park)
	mp = &MyModel_Dummy;

    /* Define here dict entries for dynamically allocated variables. */
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
    void	          *MP,	/* model parameter handle	*/
    struct tHydBrakeCU_IF *IF,	/* brake interface structure	*/
    double	           dt)	/* time step			*/
{
    int i;
    double desiredTrq[NWHEEL], dTrq[NWHEEL];
    struct tMyModel *mp = (struct tMyModel *)MP;

    for (i=0; i < NWHEEL; i++) {
	double u;

	desiredTrq[i] = mp->TrqDistrib[i] * IF->Pedal;
	dTrq[i] = desiredTrq[i] - IF->Trq_WB[i];

	if (IF->Pedal<=1e-6 && IF->Trq_WB[i]<1e-6) {
	    IF->V[i]       = 0.0;
	    IF->V[4+i]     = 0.0;
	    IF->Pedal      = 0;
	    continue;
	}

	mp->dTrqSum[i] += dTrq[i];

	/* controller P-value 1e-2; controller I-value 1e-2 */
	u = 1e-2 * dTrq[i] + 1e-2 * mp->dTrqSum[i] * dt;
	if (u>0) { 	//pressure buildup
	    if (u>1) {
		IF->V[i] = 0.0;
	    } else {
		IF->V[i] = 1.0 - u;
	    }
	    IF->V[4+i] = 0.0;
	} else {	// pressure reduction
	    IF->V[i] = 1.0;
	    if (-u>1) {
		IF->V[4+i] = 1.0;
	    } else {
		IF->V[4+i] = -u;
	    }
	}

	if (IF->V[4+i]>0)
	    IF->PumpCtrl = 1;
	else
	    IF->PumpCtrl = 0;
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
** HydBrakeCU_Register_MyModel
**
*/

int
HydBrakeCU_Register_MyModel (void)
{
    tModelClassDescr	mp;

    memset (&mp, 0, sizeof(mp));
    mp.HydBrakeControl.New = 	   	MyModel_New;
    mp.HydBrakeControl.Calc = 	   	MyModel_Calc;
    mp.HydBrakeControl.Delete =		MyModel_Delete;
    mp.HydBrakeControl.DeclQuants =	MyModel_DeclQuants;
    mp.HydBrakeControl.VersionId =	ThisVersionId;
    /* Should only be used if the model doesn't read params from extra files */
    mp.HydBrakeControl.ParamsChanged = ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_HydBrakeControl, ThisModelKind, &mp);
}
