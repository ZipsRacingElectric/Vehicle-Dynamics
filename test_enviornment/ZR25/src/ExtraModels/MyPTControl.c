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
** Simple powertrain control Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    PTControl_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "PowerTrain.Control";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


struct tMyModel {
    tGearBoxKind GBKind;
};


static void
MyModel_DeclQuants_dyn (struct tMyModel *mp, int park)
{
    static struct tMyModel MyModel_Dummy = { GBKind_NoGearBox };
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


static void
MyModel_Delete (void *MP)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Park the dict entries for dynamically allocated variables before deleting */
    MyModel_DeclQuants_dyn (mp, 1);
    free (mp);
}

/* Model output parameters in the configuration struct CfgIF, which are required
   by CarMaker, are read in before the MyModel_New() function.
   - The parametrization of these parameters is supported by the GUI.
   - These output parameters can be used internally by the model in same way like
     the input parameters
*/
static void *
MyModel_New (struct tInfos  	    *Inf,
    	     struct tPTControlCfgIF *CfgIF,
             const char     	    *KindKey)
{
    struct tMyModel *mp = NULL;
    char MsgPre[64];
    const char *ModelKind;
    int VersionId = 0;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    if (CfgIF->PTKind != PTKind_Generic) {
	LogErrF (EC_Init, "%s: supports only Generic powertrain", MsgPre);
	return NULL;
    }

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_PTControl, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*)calloc(1,sizeof(*mp));

    /* get CfgIF parameters */
    if (PTControl_GetCfgOutIF (Inf, CfgIF, ModelKind) != 0)
	return NULL;

    /* CfgIF -> Model */
    mp->GBKind = CfgIF->GearBox.GBKind;

    /* CfgIF output: verification if the parametrization corresponds to the model */
    if (CfgIF->StartEngineWithSST) {
	LogErrF (EC_Init, "%s: no support for using SST", MsgPre);
	return NULL;
    }

    return mp;
}


static int
PedalsReady2Start (struct tMyModel *mp, tPTControlIF *IF)
{
    int Ready = 0, TransmOK = 0;

    switch (mp->GBKind) {
	case (GBKind_NoGearBox):
	    TransmOK = 1;
	    break;
	case (GBKind_Manual):
	    if (IF->GearNoTrg==0 && IF->Clutch >= 0.9)
		TransmOK = 1;
	    break;
	case (GBKind_AutoWithManual):
	case (GBKind_AutoNoManual):
	    if (IF->SelectorCtrl == SelectorCtrl_N)
		TransmOK = 1;
	    break;
    }
    if (TransmOK && IF->Brake >= 0.5)
	Ready = 1;

    return Ready;
}


static int
MyModel_Calc (void *MP, tPTControlIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Operation Error */
    IF->OperationError = No_WarnError;

    /* Strategy mode */
    IF->StrategyMode = Mode_EngineDrive;

    /* Actual & Target Operation State Handling */
    switch (IF->OperationState) {
	case (OperState_Absent):
	    /* Absent -> PowerOff */
	    if (IF->Key >= KeyPos_KeyIn_PowerOff)
		IF->OperationState = OperState_PowerOff;
	    break;

	case (OperState_PowerOff):
	    /* PowerOff -> Absent */
	    if (IF->Key == KeyPos_KeyOut) {
		IF->OperationState = OperState_Absent;
		goto OutOfOperState;
	    }

	    /* PowerOff -> PowerAcc */
	    if (IF->Key >= KeyPos_KeyIn_PowerAcc) {
		IF->OperationState = OperState_PowerAccessory;
		goto OutOfOperState;
	    }
	    break;

	case (OperState_PowerAccessory):
	    /* PowerAcc -> PowerOn */
	    if (IF->Key >= KeyPos_KeyIn_PowerOn) {
		IF->Ignition       = 1;
		IF->OperationState = OperState_PowerOn;
		goto OutOfOperState;
	    }
	    break;

	case (OperState_PowerOn):
	    /* PowerOn -> PowerOff */
	    if (IF->Key<=KeyPos_KeyIn_PowerOff) {
		IF->ISGOut.Load    = 0.0;
		IF->Ignition       = 0;
		IF->OperationState = OperState_PowerOff;
		goto OutOfOperState;
	    }

	    /* PowerOn -> Driving */
	    if (PedalsReady2Start (mp, IF) && IF->Key == KeyPos_KeyIn_Starter) {
		IF->ISGOut.Load = 1.0;

		if (IF->EngineIn.Engine_on) {
		    IF->ISGOut.Load       = 0.0;
		    IF->EngineOut.set_ISC = 1;
		    IF->OperationState = OperState_Driving;
		    goto OutOfOperState;
		}
	    } else {
		IF->ISGOut.Load = 0.0;
	    }
	    break;

	case (OperState_Driving):
	    /* Driving -> PowerOn */
	    if (!IF->EngineIn.Engine_on) {
		IF->EngineOut.Load    = 0.0;
		IF->EngineOut.set_ISC = 0;
		IF->OperationState = OperState_PowerOn;
		goto OutOfOperState;
	    }

	    /* Driving -> PowerOff */
	    if (IF->Key<=KeyPos_KeyIn_PowerOff) {
		IF->EngineOut.Load    = 0.0;
		IF->EngineOut.set_ISC = 0;
		IF->Ignition          = 0;
		IF->OperationState    = OperState_PowerOff;
		goto OutOfOperState;
	    }

	    /* Gas */
	    IF->EngineOut.Load = IF->Gas;

	    break;
    }
    OutOfOperState:

    /* Clutch */
    if (mp->GBKind == GBKind_Manual)
	IF->ClutchOut.Pos = IF->Clutch;

    /* GearNo */
    if (mp->GBKind == GBKind_Manual) {
	IF->GearBoxOut.GearNoTrg = IF->GearNoTrg;
    } else if (mp->GBKind == GBKind_AutoWithManual && IF->SelectorCtrl==SelectorCtrl_M &&
	       IF->OperationState==OperState_Driving) {
	IF->GearBoxOut.GearNoTrg = IF->GearNoTrg;
    } else {
	IF->GearBoxOut.GearNoTrg = NOTSET;
    }

    return 0;
}


int
PTControl_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.PTControl.VersionId =		ThisVersionId;
    m.PTControl.New =			MyModel_New;
    m.PTControl.Calc =			MyModel_Calc;
    m.PTControl.Delete =		MyModel_Delete;
    m.PTControl.DeclQuants =		MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.PTControl.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_PTControl, ThisModelKind, &m);
}
