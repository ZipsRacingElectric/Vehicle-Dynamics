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
** Simple automatic transmission controll Model without Manumatic.
** Use this model together with the IPG "Closed" Clutch and the 
** IPG "Automatic with Converter" GearBox.
** The lock-up clutch of the converter isn't controlled in this model.
**
** For use with Examples/ExtraModels/DemoCar_MyPTTransmCU
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    TransmCU_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "PowerTrain.TCU";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


struct tMyModel {
    /* Parameters */
    double	nShiftUp_FGear[VEHICLE_NGEARS+1];
    double	nShiftDown_FGear[VEHICLE_NGEARS+1];
    double	nShiftUp_BGear[VEHICLE_NGEARS+1];
    double	nShiftDown_BGear[VEHICLE_NGEARS+1];

    /* Variables */
    int		GearNoTrg;
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


/* Model output parameters in the configuration struct CfgIF, which are required
   by CarMaker, are read in before the MyModel_New() function.
   - The parametrization of these parameters is supported by the GUI.
   - These output parameters can be used internally by the model in same way like
     the input parameters
*/
static void *
MyModel_New (struct tInfos *Inf, struct tPTTransmCU_CfgIF *CfgIF, const char *KindKey)
{
    struct tMyModel *mp = NULL;
    const char *ModelKind;
    int VersionId = 0, n, i;
    char MsgPre[64];
    double *dv = NULL;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_PTTransmCU, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*)calloc(1,sizeof(*mp));

    /* get CfgIF parameters */
    if (TransmCU_GetCfgOutIF (Inf, CfgIF, ModelKind) != 0)
	goto ErrorReturn;

    /* model will only work with a gearbox model with internal converter */
    if (CfgIF->ClKind != ClKind_Closed) {
	LogErrF (EC_Init, "%s: no 'Closed' Clutch used.", MsgPre);
	goto ErrorReturn;
    }
    if (CfgIF->GearBox.GBKind != GBKind_AutoWithManual || CfgIF->GearBox.ClKind != ClKind_Converter) {
	LogErrF (EC_Init, "%s: no 'Automatic with Converter' GearBox used.", MsgPre);
	goto ErrorReturn;
    }

    /* table to define the rot.velocity for shifting */
    n=0;
    dv = iGetTableOpt2(Inf, "MyModel.GearShift", NULL, 3, &n);

    if (dv != NULL && n > 0) {
	n = M_MIN(VEHICLE_NGEARS, n);
	for (i=0; i < n; i++) {
	    mp->nShiftDown_FGear[i+1] = dv[i + n] * rpm2radsec;
	    mp->nShiftUp_FGear[i+1]   = dv[i + 2*n] * rpm2radsec;
	    mp->nShiftDown_BGear[i+1] = dv[i + n] * rpm2radsec;
	    mp->nShiftUp_BGear[i+1]   = dv[i + 2*n] * rpm2radsec;
	}
	for (i=n; i < VEHICLE_NGEARS; i++) {
	    mp->nShiftDown_FGear[i+1] = mp->nShiftDown_FGear[n];
	    mp->nShiftUp_FGear[i+1]   = mp->nShiftUp_FGear[n];
	    mp->nShiftDown_BGear[i+1] = mp->nShiftDown_BGear[n];
	    mp->nShiftUp_BGear[i+1]   = mp->nShiftUp_BGear[n];
	}
	free(dv);
	dv = NULL;
    } else {
	LogErrF (EC_Init, "%s 'MyModel.GearShift': Wrong parametrization of GearBox-Control shift limits.", MsgPre);
	goto ErrorReturn;
    }
    mp->GearNoTrg = 0;

    /* CfgIF output: verification if the parametrization corresponds to the model */
    if (CfgIF->AutoWithMan) {
	LogErrF (EC_Init, "%s: no support for Manumatic", MsgPre);
	goto ErrorReturn;
    }

    return mp;

    ErrorReturn:
	if (mp!=NULL)
	    MyModel_Delete (mp);
	if (dv != NULL)
	    free(dv);
	return NULL;
}


static int
MyModel_Calc (void *MP, struct tPTTransmCU_IF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
#ifdef __cplusplus
    struct tPTTransmCU_IF::tPTTransmCU_IF_GearBoxOut *GBOut = &IF->GearBoxOut;
    struct tPTTransmCU_IF::tPTTransmCU_IF_GearBoxIn  *GBIn  = &IF->GearBoxIn;
#else
    struct tPTTransmCU_IF_GearBoxOut *GBOut = &IF->GearBoxOut;
    struct tPTTransmCU_IF_GearBoxIn  *GBIn  = &IF->GearBoxIn;
#endif
    const tPTTransmCU_CfgIF *CfgIF  = IF->CfgIF;

    /* if the powertrain is not running, do nothing */
    if (!IF->Ignition)
	goto ReturnZero;

    /* always open, the Lock-up function of the Converter is not used */
    GBOut->Clutch.Pos = 1.0;

    /* use the GearNoTrg from PTControl, if it is set */
    if (GBOut->GearNoTrg != NOTSET)
	mp->GearNoTrg = M_BOUND(-CfgIF->GearBox.nBGears, CfgIF->GearBox.nFGears, GBOut->GearNoTrg);

    /* calculate the GearNoTrg */
    if (IF->SelectorCtrl != SelectorCtrl_M && GBOut->GearNoTrg == NOTSET) {
	double RotvRef, nUp=0, nDown=0;
	
	/* get the input rot.velocity into the GearBox as a reference velocity */
	RotvRef = fabs(GBIn->Clutch.rotv_out);
	
	/* getting rot.velocity when to shift depending on the gear number */
	nUp   = mp->nShiftUp_FGear[abs(GBIn->GearNo)];
	nDown = mp->nShiftDown_FGear[abs(GBIn->GearNo)];
	
	if (IF->SelectorCtrl == SelectorCtrl_D) {	/* driving forwards */
	    if (GBIn->GearNo > 0) {	/* already driving forwards */
		/* upshifting */
		if((GBIn->GearNo < CfgIF->GearBox.nFGears) && (RotvRef >= nUp))
		    mp->GearNoTrg = GBIn->GearNo + 1;
		/* downshifting */
		if((GBIn->GearNo > 1) && (RotvRef <= nDown))
		    mp->GearNoTrg = GBIn->GearNo - 1;
	    } else {			/* was driving backwards or in neutral position */
		mp->GearNoTrg = 1;
	    }
	} else if (IF->SelectorCtrl == SelectorCtrl_R) {	/* driving backwards */
	    if (GBIn->GearNo < 0) {	/* already driving backwards */
		if((-GBIn->GearNo < CfgIF->GearBox.nBGears) && (RotvRef >= nUp))
		    mp->GearNoTrg = GBIn->GearNo - 1;
		if((-GBIn->GearNo > 1) && (RotvRef <= nDown))
		    mp->GearNoTrg = GBIn->GearNo + 1;
	    } else {			/* was driving forwards or in neutral position */
		mp->GearNoTrg = -1;
	    }
	} else if (IF->SelectorCtrl == SelectorCtrl_P) {	/* parking position */
	    mp->GearNoTrg     = 0;
	    GBOut->set_ParkBrake = 1;
	} else {
	    mp->GearNoTrg     = 0;
	}
    }

    GBOut->GearNoTrg = mp->GearNoTrg;

    return 0;

    ReturnZero:
	GBOut->Clutch.Pos    = 1.0;
	GBOut->GearNoTrg     = 0;
	GBOut->set_ParkBrake = 0;
	return 0;
}


int 
TransmCU_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.PTTransmCU.VersionId =	ThisVersionId;
    m.PTTransmCU.New =		MyModel_New;
    m.PTTransmCU.Calc =		MyModel_Calc;
    m.PTTransmCU.Delete =	MyModel_Delete;
    m.PTTransmCU.DeclQuants =	MyModel_DeclQuants;
    m.PTTransmCU.PreSimSetup =	NULL;
    /* Should only be used if the model doesn't read params from extra files */
    m.PTTransmCU.ParamsChanged = ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_PTTransmCU, ThisModelKind, &m);
}
