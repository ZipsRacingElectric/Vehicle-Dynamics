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
** Simple engine controll Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    EngineCU_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


struct tMyModel {
    struct {
    	double I;
    	double P;
    	double Load_I;
    	double Load_P;
    } ISCtrl;

    double	nFuelCutOff;
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



static void *
MyModel_New (struct tInfos *Inf, struct tPTEngineCU_CfgIF *CfgIF, const char *KindKey)
{
    struct tMyModel *mp = NULL;
    const char *ModelKind;
    int VersionId = 0;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_PTEngineCU, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*)calloc(1,sizeof(*mp));

    /* Idle Speed Controller */
    mp->ISCtrl.P = iGetDblOpt (Inf, "PowerTrain.MyEngineCU.ISCtrl.P", 	  0.1);
    mp->ISCtrl.I = iGetDblOpt (Inf, "PowerTrain.MyEngineCU.ISCtrl.I", 	  50);

    /* FuelCutOff */
    mp->nFuelCutOff = iGetDblOpt(Inf, "PowerTrain.ECU.FuelCutOff",
					CfgIF->rotv_idle * 2.0 * radsec2rpm)*rpm2radsec;

    return mp;
}


static int
MyModel_Calc (void *MP, struct tPTEngineCU_IF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    const tPTEngineCU_CfgIF *CfgIF = IF->CfgIF;

    if (!IF->Ignition)
	goto ReturnZero;

    /* Engine on ? */
    if (IF->rotv > CfgIF->rotv_idle)
	IF->Engine_on = 1;
    if (IF->rotv < CfgIF->rotv_off)
	IF->Engine_on = 0;

    if (IF->Engine_on) {
	/** TrqFull & TrqDrag & TrqOpt */
	if (CfgIF->TrqFull != NULL)
	    IF->TrqFull = LMEval(CfgIF->TrqFull, IF->rotv);
	if (CfgIF->TrqDrag != NULL)
	    IF->TrqDrag = LMEval(CfgIF->TrqDrag, IF->rotv);
	if (CfgIF->TrqOpt != NULL)
	    IF->TrqOpt = LMEval(CfgIF->TrqOpt, IF->rotv);

	/** Engine Load */
	if (IF->Load == NOTSET)
	    IF->Load = 0.0;
	
	/** Idle Speed Control */
	if (IF->set_ISC) {
	    if (IF->rotv < CfgIF->rotv_idle + 200.0*rpm2radsec) {
		double Load_I, drotv, Load;
		drotv = CfgIF->rotv_idle - IF->rotv;

		if (IF->Load < 1e-3) {
		    mp->ISCtrl.Load_I += drotv * mp->ISCtrl.I * dt;
		    mp->ISCtrl.Load_I  = M_BOUND (0.0, 1.0, mp->ISCtrl.Load_I);
		    Load_I = mp->ISCtrl.Load_I;
		} else {
		    double Ifac;
		    Ifac = 1.0 + drotv / (200.0*rpm2radsec);
		    Ifac = M_BOUND (0.0, 1.0, Ifac);
		    Load_I = mp->ISCtrl.Load_I * Ifac;
		}

		mp->ISCtrl.Load_P = M_MAX(0.0, drotv * mp->ISCtrl.P);

		Load = Load_I + mp->ISCtrl.Load_P;
		Load = M_BOUND (0.0, 1.0, Load);
		IF->Load += Load;
	    } else {
		mp->ISCtrl.Load_I = mp->ISCtrl.Load_P = 0.0;
	    }
	}	
	IF->Load = M_BOUND (0.0, 1.0, IF->Load);

	/** FuelCutOff */
	if (IF->FuelCutOff == NOTSET) {
	    if (IF->rotv >= mp->nFuelCutOff && IF->Load<=1e-3)
		IF->FuelCutOff = 1;
	    else
		IF->FuelCutOff = 0;

	} else if (IF->FuelCutOff == 1) {
	    IF->Load            = 0.0;
	}
    } else {
	goto ReturnZero;
    }


    return 0;

    ReturnZero:
	IF->Engine_on  = 0;
	IF->FuelCutOff = 0;
	IF->Load       = 0.0;
	IF->TrqDrag    = 0.0;
	IF->TrqFull    = 0.0;
	IF->TrqOpt     = 0.0;
	return 0;
}


int 
EngineCU_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.PTEngineCU.VersionId =		ThisVersionId;
    m.PTEngineCU.New =			MyModel_New;
    m.PTEngineCU.Calc =			MyModel_Calc;
    m.PTEngineCU.Delete =		MyModel_Delete;
    m.PTEngineCU.DeclQuants =		MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.PTEngineCU.ParamsChanged =	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_PTEngineCU, ThisModelKind, &m);
}
