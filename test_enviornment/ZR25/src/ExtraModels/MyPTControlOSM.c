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
** Simple operation state machine Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    PTControlOSM_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "PowerTrain.ControlOSM";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


struct tMyModel {
    tGearBoxKind GBKind;
    int 	 State;
    int 	 Func_used[OSMFunc_MaxNumber];
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

    if (mp!=NULL)
	free (mp);
    mp = NULL;
}


static void *
MyModel_New (struct tInfos  	    	*Inf,
    	     struct tPTControlOSM_CfgIF *CfgIF,
             const char     	    	*KindKey)
{
    struct tMyModel *mp = NULL;
    char MsgPre[64];
    const char *ModelKind;
    int VersionId = 0;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_PTControlOSM, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    if (CfgIF->StartEngineWithSST) {
	LogErrF (EC_Init, "%s: supports only starting with key", MsgPre);
	return NULL;
    }
    mp = (struct tMyModel*)calloc(1,sizeof(*mp));

    /* CfgIF -> Model */
    mp->GBKind = CfgIF->GBKind;
    memcpy (mp->Func_used, CfgIF->Func_used, OSMFunc_MaxNumber*sizeof(int));

    if (!mp->Func_used[OSMFunc_Drive]) {
	LogErrF (EC_Init, "%s: no function for drive exists", MsgPre);
	goto ErrorReturn;
    }

    return mp;

    ErrorReturn:
	MyModel_Delete (mp);
	return NULL;
}


static int
PedalsReady2Start (struct tMyModel *mp, tPTControlOSM_IF *IF)
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
MyModel_Calc_Absent (struct tMyModel *mp, tPTControlOSM_IF *IF, double dt)
{
    int rv = 0;

    switch (mp->State) {
	case 0: /* Absent Function exist ? */
	    if (mp->Func_used[OSMFunc_Absent]) {
		IF->OSMFunc = OSMFunc_Absent;
		mp->State++;
	    } else {
		mp->State+=2;
	    }
	    break;

	case 1: /* Wait for Absent is ready */
	    if (IF->FuncReturn)
		mp->State++;
	    break;

	case 2: /* Absent -> PowerOff ? */
	    IF->OSMFunc = OSMFunc_None;
	    if (IF->Key >= KeyPos_KeyIn_PowerOff)
		mp->State++;
	    break;

	case 3: /* PowerOff function exist ? */
	    if (mp->Func_used[OSMFunc_PowerOff]) {
		IF->OSMFunc = OSMFunc_PowerOff;
		mp->State++;
	    } else {
		mp->State+=2;
	    }
	    break;

	case 4: /* Wait for PowerOff is ready */
	    if (IF->FuncReturn)
		mp->State++;
	    break;

	case 5: /* PowerOff reached */
	    IF->OperationState = OperState_PowerOff;
	    rv = 1;
	    break;
    }

    return rv;
}


static int
MyModel_Calc_PowerOff (struct tMyModel *mp, tPTControlOSM_IF *IF, double dt)
{
    int rv = 0;

     /* PowerOff -> Absent ? */
    if (IF->Key == KeyPos_KeyOut) {
	IF->OperationState = OperState_Absent;
	return rv=1;
    }

    switch (mp->State) {
	case 0:	/* PowerOff -> PowerAcc ? */
	    if (IF->Key >= KeyPos_KeyIn_PowerAcc)
		mp->State++;
	    break;
	
	case 1: /* PowerAcc function exist ? */
	    if (mp->Func_used[OSMFunc_PowerAcc]) {
		IF->OSMFunc = OSMFunc_PowerAcc;
		mp->State++;
	    } else {
		mp->State+=2;
	    }
	    break;

	case 2: /* Wait for PowerAcc is ready */
	    if (IF->FuncReturn)
		mp->State++;
	    break;

	case 3: /* PowerAcc reached */
	    IF->OperationState = OperState_PowerAccessory;
	    rv = 1;
	    break;
    }

    return rv;
}


static int
MyModel_Calc_PowerAcc (struct tMyModel *mp, tPTControlOSM_IF *IF, double dt)
{
    int rv = 0;

    switch (mp->State) {
	case 0:	/* PowerAcc -> PowerOn ? */
	    if (IF->Key >= KeyPos_KeyIn_PowerOn)
		mp->State++;
	    break;
	
	case 1: /* PowerOn function exist ? */
	    if (mp->Func_used[OSMFunc_PowerOn]) {
		IF->OSMFunc = OSMFunc_PowerOn;
		mp->State++;
	    } else {
		mp->State+=2;
	    }
	    break;

	case 2: /* Wait for PowerOn is ready */
	    if (IF->FuncReturn)
		mp->State++;
	    break;

	case 3: /* PowerOn reached */
	    IF->OperationState = OperState_PowerOn;
	    rv = 1;
	    break;
    }

    return rv;
}


static int
MyModel_Calc_PowerOn (struct tMyModel *mp, tPTControlOSM_IF *IF, double dt)
{
    int rv = 0;

    switch (mp->State) {
	case 0: /* PowerOn -> Driving or PowerOn -> PowerOff */

	    if (PedalsReady2Start (mp, IF)) {
		/* PowerOn -> Driving */
		if (IF->Key == KeyPos_KeyIn_Starter) {
		    mp->State++;
		}
	    } else {
		/* PowerOn -> PowerOff */
		if (IF->Key<=KeyPos_KeyIn_PowerOff)
		    mp->State  = 9;
	    }
	    break;

	case 1: /* Start function exist ? */
	    if (mp->Func_used[OSMFunc_Start]) {
		IF->OSMFunc = OSMFunc_Start;
		mp->State++;
	    } else {
		mp->State+=2;
	    }
	    break;

	case 2: /* Wait for Start is ready */
	    /* Skip Start if Key is not Starter */
	    if (IF->Key != KeyPos_KeyIn_Starter)
		mp->State = 4;

	    if (IF->FuncReturn)
		mp->State++;
	    break;

	case 3: /* Drive reached */
	    IF->OperationState = OperState_Driving;
	    rv = 1;
	    break;

	case 4: /* SkipStart function exist ? */
	    if (mp->Func_used[OSMFunc_SkipStart]) {
		IF->OSMFunc = OSMFunc_SkipStart;
		mp->State++;
	    } else {
		mp->State+=2;
	    }
	    break;

	case 5: /* Wait for SkipStart is ready */
	    if (IF->FuncReturn)
		mp->State++;
	    break;

	case 6: /* SkipStart reached -> PowerOn function exist ? */
	    if (mp->Func_used[OSMFunc_PowerOn]) {
		IF->OSMFunc = OSMFunc_PowerOn;
		mp->State++;
	    } else {
		mp->State+=2;
	    }
	    break;

	case 7: /* Wait for PowerOn is ready */
	    if (IF->FuncReturn)
		mp->State++;
	    break;

	case 8: /* PowerOn reached -> back to beginning */
	    IF->OperationState = OperState_PowerOn;
	    rv = 1;
	    break;

	case 9: /* PowerOff function exist ? */
	    if (mp->Func_used[OSMFunc_PowerOff]) {
		IF->OSMFunc = OSMFunc_PowerOff;
		mp->State++;
	    } else {
		mp->State+=2;
	    }
	    break;

	case 10: /* Wait for PowerOff is ready */
	    if (IF->FuncReturn)
		mp->State++;
	    break;

	case 11: /* PowerOff reached */
	    IF->OperationState = OperState_PowerOff;
	    rv = 1;
	    break;
    }

    return rv;
}


static int
MyModel_Calc_Drive (struct tMyModel *mp, tPTControlOSM_IF *IF, double dt)
{
    int rv = 0;

    switch (mp->State) {
	case 0: /* Set Drive */
	    IF->OSMFunc = OSMFunc_Drive;
	    mp->State++;
	    break;

	case 1: /* Treat error or Drive -> PowerOff */
	    if (IF->FuncReturn < 0) {
		mp->State++;
		break;
	    }

	    /* Drive -> PowerOff */
	    if (IF->Key<=KeyPos_KeyIn_PowerOff)
		mp->State = 5;
	    break;

	case 2: /* PowerOn function exist ? */
	    if (mp->Func_used[OSMFunc_PowerOn]) {
		IF->OSMFunc = OSMFunc_PowerOn;
		mp->State++;
	    } else {
		mp->State+=2;
	    }
	    break;

	case 3: /* Wait for PowerOn is ready */
	    if (IF->FuncReturn)
		mp->State++;
	    break;

	case 4: /* PowerOn reached */
	    IF->OperationState = OperState_PowerOn;
	    rv = 1;
	    break;

	case 5: /* PowerOff function exist ? */
	    if (mp->Func_used[OSMFunc_PowerOff]) {
		IF->OSMFunc = OSMFunc_PowerOff;
		mp->State++;
	    } else {
		mp->State+=2;
	    }
	    break;

	case 6: /* Wait for PowerOff is ready */
	    if (IF->FuncReturn)
		mp->State++;
	    break;

	case 7: /* PowerOff reached */
	    IF->OperationState = OperState_PowerOff;
	    rv = 1;
	    break;
    }

    return rv;
}


static int
MyModel_Calc (void *MP, tPTControlOSM_IF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    int rv = 0;

    switch (IF->OperationState) {
	case (OperState_Absent):
	    rv = MyModel_Calc_Absent (mp, IF, dt);	break;

	case (OperState_PowerOff):
	    rv = MyModel_Calc_PowerOff (mp, IF, dt);	break;

	case (OperState_PowerAccessory):
	    rv = MyModel_Calc_PowerAcc (mp, IF, dt);	break;

	case (OperState_PowerOn):
	    rv = MyModel_Calc_PowerOn (mp, IF, dt);	break;

	case (OperState_Driving):
	    rv = MyModel_Calc_Drive (mp, IF, dt);	break;
    }

    /* Operation reached ? */
    if (rv == 1) {
	IF->OSMFunc = OSMFunc_None;
	mp->State   = 0;
    }

    return 0;
}


int
PTControlOSM_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.PTControlOSM.VersionId =		ThisVersionId;
    m.PTControlOSM.New =		MyModel_New;
    m.PTControlOSM.Calc =		MyModel_Calc;
    m.PTControlOSM.Delete =		MyModel_Delete;
    m.PTControlOSM.DeclQuants =		MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.PTControlOSM.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_PTControlOSM, ThisModelKind, &m);
}
