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
** Simple Vehicle Operator Model to start powertrain with ignition key
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    VhclOperator_Register_MyModel ()
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

typedef struct {
    tKeyPosition Key;
    int	GearNo;
    int	SelectorCtrl;
    double	Gas;
    double	Clutch;
    double	Brake;
    double	BrakePark;
} tIPGOperator_State;

/* Model Parameters */
struct tMyModel {
    tGearBoxKind	GearBoxKind;
    int			Phase;
    tOperationState	OperationState_last;
    tOperationState	OperationState_trg_last;
    tIPGOperator_State	Absent, PowerOff, PowerAcc, PowerOn, Starting, Driving, Braking;
};

static const char ThisModelClass[] = "DrivMan.VhclOperator";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


/******************************************************************************/


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


static int
MyModel_PreSimSetup (void *MP, struct tVhclOperatorPreSimIF *PreSimIF,
		  	       struct tVhclOperatorIF       *IF)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Take testrun start conditions only for state Driving */
    if (IF->OperationState != OperState_Driving)
	return 0;

    if (mp->GearBoxKind==GBKind_Manual) {
	IF->GearNo = PreSimIF->GearNo;
	IF->Clutch = PreSimIF->Clutch;
    } else {
	IF->SelectorCtrl = PreSimIF->SelectorCtrl;
    }
    IF->Steering.Ang = PreSimIF->SteerAng;

    return 0;
}


static void
SetOutput (struct tMyModel *mp, tIPGOperator_State *State, tVhclOperatorIF *IF)
{
    if (mp->GearBoxKind==GBKind_Manual) {
	IF->GearNo = State->GearNo;
	IF->Clutch = State->Clutch;
    } else {
	IF->SelectorCtrl = State->SelectorCtrl;
    }
    IF->Key        = State->Key;
    IF->Gas        = State->Gas;
    IF->Brake      = State->Brake;
    IF->BrakePark  = State->BrakePark;
}

static void
SetOutputAndWait (struct tMyModel *mp, tIPGOperator_State *State, tVhclOperatorIF *IF)
{
    switch (mp->Phase) {
	case 0:
	    SetOutput (mp, State, IF);
	    mp->Phase++;
	    break;
	case 1:
	    /* Wait for target Operation State */
	    break;
    }
}


static void
ChangeState_Drive2Off (struct tMyModel *mp, tVhclOperatorIF *IF)
{
    switch (mp->Phase) {
	case 0:
	    /* Brake to standstill ? */
	    if (IF->Velocity > 0.1)
		mp->Phase++;
	    else
		mp->Phase = 3;
	    break;
	case 1:
	    /* Active brake */
	    SetOutput (mp, &mp->Braking, IF);
	    mp->Phase++;
	    break;
	case 2:
	    /* Wait until standstill */
	    if (IF->Velocity <= 0.1)
		mp->Phase++;
	    break;
	case 3:
	    IF->Key          = mp->PowerOff.Key;
	    IF->Steering.Ang = 0.0;
	    mp->Phase++;
	    break;
	case 4:
	    /* Wait for PowerOff */
	    break;
    }
}


static int
MyModel_Calc (void *MP, tVhclOperatorIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    if (IF->OperationState     != mp->OperationState_last
     || IF->OperationState_trg != mp->OperationState_trg_last) {
	mp->Phase = 0;
    }

    switch (IF->OperationState_trg) {
	/** Target state: Absent */
	case (OperState_Absent):
	    switch (IF->OperationState) {
		/* Current state: Absent */
		case (OperState_Absent):
		    SetOutput (mp, &mp->Absent, IF);
		    goto Finished;

		/* Current state: PowerOff */
		case (OperState_PowerOff):
		    SetOutputAndWait (mp, &mp->Absent, IF);
		    break;
		
		/* Current state: PowerAccessory */
		case (OperState_PowerAccessory):
		    SetOutputAndWait (mp, &mp->PowerOn, IF);
		    break;

		/* Current state: PowerOn */
		case (OperState_PowerOn):
		    SetOutputAndWait (mp, &mp->PowerOff, IF);
		    break;

		/* Current state: Driving */
		case (OperState_Driving):
		    ChangeState_Drive2Off (mp, IF);
		    break;
	    }
	    break;

	/** Target state: PowerOff */
	case (OperState_PowerOff):
	    switch (IF->OperationState) {
		/* Current state: Absent */
		case (OperState_Absent):
		    SetOutput (mp, &mp->PowerOff, IF);
		    break;

		/* Current state: PowerOff */
		case (OperState_PowerOff):
		    SetOutput (mp, &mp->PowerOff, IF);
		    goto Finished;
		
		/* Current state: PowerAccessory */
		case (OperState_PowerAccessory):
		    SetOutputAndWait (mp, &mp->PowerOn, IF);
		    break;

		/* Current state: PowerOn */
		case (OperState_PowerOn):
		    SetOutputAndWait (mp, &mp->PowerOff, IF);
		    break;

		/* Current state: Driving */
		case (OperState_Driving):
		    ChangeState_Drive2Off (mp, IF);
		    break;
	    }
	    break;

	/** Target state: PowerAccessory */
	case (OperState_PowerAccessory):
	    switch (IF->OperationState) {
		/* Current state: Absent */
		case (OperState_Absent):
		    SetOutput (mp, &mp->PowerOff, IF);
		    break;

		/* Current state: PowerOff */
		case (OperState_PowerOff):
		    SetOutputAndWait (mp, &mp->PowerAcc, IF);
		    break;
		
		/* Current state: PowerAccessory */
		case (OperState_PowerAccessory):
		    SetOutput (mp, &mp->PowerAcc, IF);
		    goto Finished;

		/* Current state: PowerOn */
		case (OperState_PowerOn):
		    SetOutputAndWait (mp, &mp->PowerOff, IF);
		    break;

		/* Current state: Driving */
		case (OperState_Driving):
		    ChangeState_Drive2Off (mp, IF);
		    break;
	    }
	    break;

	/** Target state: PowerOn */
	case (OperState_PowerOn):
	    switch (IF->OperationState) {
		/* Current state: Absent */
		case (OperState_Absent):
		    SetOutput (mp, &mp->PowerOff, IF);
		    break;

		/* Current state: PowerOff */
		case (OperState_PowerOff):
		    SetOutputAndWait (mp, &mp->PowerAcc, IF);
		    break;
		
		/* Current state: PowerAccessory */
		case (OperState_PowerAccessory):
		    SetOutputAndWait (mp, &mp->PowerOn, IF);
		    break;

		/* Current state: PowerOn */
		case (OperState_PowerOn):
		    SetOutput (mp, &mp->PowerOn, IF);
		    goto Finished;

		/* Current state: Driving */
		case (OperState_Driving):
		    ChangeState_Drive2Off (mp, IF);
		    break;
	    }
	    break;

	/** Target state: Driving */
	case (OperState_Driving):
	    switch (IF->OperationState) {
		/* Current state: Absent */
		case (OperState_Absent):
		    SetOutput (mp, &mp->PowerOff, IF);
		    break;

		/* Current state: PowerOff */
		case (OperState_PowerOff):
		    SetOutputAndWait (mp, &mp->PowerAcc, IF);
		    break;
		
		/* Current state: PowerAccessory */
		case (OperState_PowerAccessory):
		    SetOutputAndWait (mp, &mp->PowerOn, IF);
		    break;

		/* Current state: PowerOn */
		case (OperState_PowerOn):
		    SetOutputAndWait (mp, &mp->Starting, IF);
		    break;

		/* Current state: Driving */
		case (OperState_Driving):
		    SetOutput (mp, &mp->Driving, IF);
		    goto Finished;
	    }
	    break;
    }

    IF->OperatorFinished        = 0;
    mp->OperationState_last     = IF->OperationState;
    mp->OperationState_trg_last = IF->OperationState_trg;
    return 0;

    Finished:
	IF->OperatorFinished        = 1;
	mp->OperationState_last     = IF->OperationState;
	mp->OperationState_trg_last = IF->OperationState_trg;
	mp->Phase                   = 0;
	return 0;
}


static void *
MyModel_New (
    struct tInfos 		*Inf,
    struct tVhclOperatorCfgIF   *CfgIF,
    const char 			*KindKey)
{
    struct tMyModel *mp = NULL;
    tIPGOperator_State *State = NULL;
    char MsgPre[64];
    const char *ModelKind;
    int VersionId = 0;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_VhclOperator, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    if (CfgIF->StartEngineWithSST) {
	LogErrF(EC_Init, "%s: starting with start-stop button not supported", MsgPre);
	return NULL;
    }

    mp = (struct tMyModel *)calloc(1, sizeof(*mp));

    /* CfgIF -> Model */
    mp->GearBoxKind = CfgIF->GearBoxKind;

    /* Model Parameters */
    State = &mp->Absent;
    State->Key          = KeyPos_KeyOut;
    State->SelectorCtrl = SelectorCtrl_P;
    State->BrakePark    = 1.0;

    State = &mp->PowerOff;
    State->Key          = KeyPos_KeyIn_PowerOff;
    State->SelectorCtrl = SelectorCtrl_P;
    State->BrakePark    = 1.0;

    State = &mp->PowerAcc;
    State->Key          = KeyPos_KeyIn_PowerAcc;
    State->SelectorCtrl = SelectorCtrl_P;
    State->BrakePark    = 1.0;

    State = &mp->PowerOn;
    State->Key          = KeyPos_KeyIn_PowerOn;
    State->SelectorCtrl = SelectorCtrl_P;
    State->BrakePark    = 1.0;

    State = &mp->Starting;
    State->Key          = KeyPos_KeyIn_Starter;
    State->SelectorCtrl = SelectorCtrl_N;
    State->Clutch       = 1.0;
    State->Brake        = 1.0;

    State = &mp->Driving;
    State->Key          = KeyPos_KeyIn_PowerOn;
    State->GearNo       = 1;
    State->SelectorCtrl = SelectorCtrl_D;
    State->Clutch       = 1.0;

    State = &mp->Braking;
    State->Key          = KeyPos_KeyIn_PowerOn;
    State->SelectorCtrl = SelectorCtrl_N;
    State->Clutch       = 1.0;
    State->Brake        = 0.6;

    return mp;
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
VhclOperator_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.VhclOperator.VersionId =	 	ThisVersionId;
    m.VhclOperator.New    = 		MyModel_New;
    m.VhclOperator.Calc   = 		MyModel_Calc;
    m.VhclOperator.DeclQuants   = 	MyModel_DeclQuants;
    m.VhclOperator.Delete = 		MyModel_Delete;
    m.VhclOperator.PreSimSetup = 	MyModel_PreSimSetup;
    /* Should only be used if the model doesn't read params from extra files */
    m.VhclOperator.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_VhclOperator, ThisModelKind, &m);
}
