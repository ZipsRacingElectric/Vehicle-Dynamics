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
** Simple powertrain Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    PowerTrain_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "Car/PowerFlow.h"
#include "MyModels.h"

#define NWHEEL 4

static const char ThisModelClass[] = "PowerTrain";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


struct tMyModel {
    /* Parameters */
    double Gas2Trq;	/* Coefficient Gas -> DriveTorque at Wheel */

    struct tMyWheel {
	double	Irot;
	double	Irot_act;
	double	rota;
    } Whl[NWHEEL];
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

/* Model output parameters in the configuration struct CfgIF, which are required
   by CarMaker, are read in before the MyModel_New() function.
   - The parametrization of these parameters is supported by the GUI.
   - These output parameters can be used internally by the model in same way like
     the input parameters
*/
static void *
MyModel_New (tInfos *Inf, struct tPowerTrainCfgIF *CfgIF, const char *KindKey)
{
    struct tMyModel *mp = NULL;
    int i, VersionId = 0;
    char MsgPre[64];
    const char *ModelKind;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_PowerTrain, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*)calloc(1,sizeof(*mp));

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    /* get CfgIF parameters */
    if (PowerTrain_GetCfgOutIF (Inf, CfgIF, ModelKind) != 0)
	goto ErrorReturn;

    if (CfgIF->nWheels != NWHEEL) {
	LogErrF(EC_Init, "%s: model supports only a four wheel vehicle", MsgPre);
	goto ErrorReturn;
    }

    /* Gas -> Drive Torque */
    mp->Gas2Trq = iGetDbl(Inf, "MyPowerTrain.Gas2Trq");

    for (i=0; i<NWHEEL; i++)
	mp->Whl[i].Irot = CfgIF->Wheel_Iyy[i];

    /* CfgIF output: verification if the parametrizationcorresponds to the model */
    if (CfgIF->PTKind != PTKind_BEV) {
	LogErrF (EC_Init, "%s: model supports only electrical powertrain", MsgPre);
	goto ErrorReturn;
    }
    if (CfgIF->GBKind != GBKind_NoGearBox) {
	LogErrF (EC_Init, "%s: model supports only no gearbox", MsgPre);
	goto ErrorReturn;
    }
    if (CfgIF->StartEngineWithSST) {
	LogErrF (EC_Init, "%s: no support for starting with SST", MsgPre);
	goto ErrorReturn;
    }

    return mp;

    ErrorReturn:
	free (mp);
	return NULL;
}


static int
MyModel_PreSimSetup (void *MP, struct tPowerTrainPreSimIF *PreSimIF,
			       struct tPowerTrainIF       *IF)
{
    int iS;

    /* Initial wheel speed */
    for (iS=0; iS<NWHEEL; iS++)
	IF->WheelOut[iS].rotv = PreSimIF->Whl_rotv[iS];

    return 0;
}


static int
MyModel_Calc (void *MP, struct tPowerTrainIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    int i;
    double Irot[NWHEEL], Irot_act[NWHEEL], Trq_Ext2W[NWHEEL];

    /* PowerTrain Operation State */
    if (SimCore.State == SCState_Simulate) {
        /* In simulation phase powertrain is always on/ready */
	IF->OperationState = OperState_Driving;
    } else {
        /* In non simulation phase use substitute model:
	   Operation State corresponds to target Operation State */
	IF->OperationState = PowerTrain_TargetOperationState ();
    }

    /* Calculate the drive torque at wheels */
    for (i=0; i<NWHEEL; i++)
	IF->WheelOut[i].Trq_Drive = IF->Gas * mp->Gas2Trq;

    /* Support torque to the wheel carrier */
    for (i=0; i<NWHEEL; i++)
	IF->WheelOut[i].Trq_Supp2WC = -IF->WheelOut[i].Trq_Drive;

    /* Reduce the brake torques and calculate the wheel inertia for low
       velocity */
    for (i=0; i<NWHEEL; i++)
	Irot[i] = mp->Whl[i].Irot;
    PT_Wheels_UpdateTrq (dt, Irot, Irot_act, Trq_Ext2W);
    for (i=0; i<NWHEEL; i++)
	mp->Whl[i].Irot_act = Irot_act[i];

    /* Wheel accelerations */
    for (i=0; i<NWHEEL; i++)
	mp->Whl[i].rota = (IF->WheelOut[i].Trq_Drive + Trq_Ext2W[i]) /
	    mp->Whl[i].Irot_act;

    /* Integration */
    if (SimCore.State != SCState_Simulate && SimCore.State != SCState_EndIdleGet) {
	for (i=0; i<NWHEEL; i++)
	    mp->Whl[i].rota = 0.0;
    }
    if (SimCore.State == SCState_EndIdleSet) {
	for (i=0; i<NWHEEL; i++)
	    IF->WheelOut[i].rotv = 0.0;
    } else {
	for (i=0; i<NWHEEL; i++) {
	    IF->WheelOut[i].rotv += mp->Whl[i].rota      * dt;
	    IF->WheelOut[i].rot  += IF->WheelOut[i].rotv * dt;
	}
    }

    /* Support torque to the vehicle body */
    IF->Trq_Supp2Bdy1[1]  = 0.0;
    IF->Trq_Supp2Bdy1B[1] = 0.0;

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
PowerTrain_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.PowerTrain.VersionId =		ThisVersionId;
    m.PowerTrain.New =			MyModel_New;
    m.PowerTrain.Calc =			MyModel_Calc;
    m.PowerTrain.Delete =		MyModel_Delete;
    m.PowerTrain.DeclQuants =		MyModel_DeclQuants;
    m.PowerTrain.PreSimSetup =		MyModel_PreSimSetup;
    /* Should only be used if the model doesn't read params from extra files */
    m.PowerTrain.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_PowerTrain, ThisModelKind, &m);
}
