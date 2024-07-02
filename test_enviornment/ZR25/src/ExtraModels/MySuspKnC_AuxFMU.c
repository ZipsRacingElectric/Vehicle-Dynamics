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
** Auxiliary model for "MySuspKnC_FMU"
**
** For every suspension kinematics FMU an auxiliary model has to be provided
** - the examples from "MySuspKnC.c" are used as auxiliary models here
** - model is used whenever SimCore.State!=SCState_Simulate, e.g. for the
**     calculation of start conditions and model check
** - the auxiliary model behaviour has to be similar to the FMU model
**     behaviour under static conditions
** - model has to be capable of ForceCoupling if SuspKnC.ForceCoupling==1
** - model name has to start with an underscore
** - provide your own matching auxiliary models if
** 	- you use your own FMU or
** 	- you changed parameters in the example FMU
**
******************************************************************************
*/

#include "Global.h"
#include "Log.h"
#include "InfoUtils.h"
#include "SimCore.h"
#include "FMUAuxModel.h"

#include "Vehicle/MBSUtils.h"
#include "Car/Susp.h"

#include "MySuspKnC.c"

static const char ThisModelKindAuxFMU[]   = "_MyModel_AuxFMU";

/* Model parameter */
typedef struct tAuxData{
    /* Mandatory first member, contains extra parameters for aux_Calc(). */
    tFMU_AuxHead	Head;		/* Do not change */

    /* Pointer to auxiliary model
     * - Replace "tMyModel" with user specific struct
     * - do not change the name "MyAuxModel" */
    struct tMyModel	*MyAuxModel;
} tAuxData;

static void *
MyModel_New_AuxFMU (
	tInfos  	*Inf,
	tSuspCfgIF	*SuspCfgIF,
	const char     	*KindKey,
	const char     	*Pre)
{
    struct tAuxData *mp = NULL;
    mp = (struct tAuxData*)calloc(1, sizeof(*mp));

    /* Replace this block with user model if FMU other than MySuspKnC_FMU is used
     * - replace "tMyModel" with user specific struct
     * - do not change the name "MyAuxModel" */
    {
	if (SuspCfgIF->ForceCoupling) {
	    mp->MyAuxModel = (struct tMyModel *) MyModel_New_FrcCpl(Inf, SuspCfgIF, KindKey, Pre);
	} else {
	    mp->MyAuxModel = (struct tMyModel *) MyModel_New(Inf, SuspCfgIF, KindKey, Pre);
	}
    }

    return mp;
}

static int
MyModel_Calc_AuxFMU(void *MP, tSuspIF *IFMain, tSuspIF *IFOpp, double dt)
{
    tAuxData *mp = (struct tAuxData*)MP;

    if (mp->Head.CallType == 0) {
	LogErrF(EC_Sim, "Model '%s' must not be used outside the context of an FMU", ThisModelKindAuxFMU);
	return -1;
    }

    if (SimCore.State == SCState_Simulate)
	return FMU_AuxAction_Pre_FMU;	/* Invoke FMU afterwards!" */

    /* Replace this block with user model if FMU other than MySuspKnC_FMU is used */
    {
	if (IFMain->CfgIF->ForceCoupling) {
	    MyModel_Calc_FrcCpl(mp->MyAuxModel, IFMain, IFOpp, dt);
	} else {
	    MyModel_Calc(mp->MyAuxModel, IFMain, IFOpp, dt);
	    double l0 = 0.24;           /* Length offset ForceElements */
	    IFMain->Kin[ixSpring] += l0;
	    IFMain->Kin[ixDamp]   += l0;
	    IFMain->Kin[ixBuf]    += l0;
	}
    }

    return FMU_AuxAction_Pre;		/* "Don't invoke FMU afterwards!" */
}

static void
MyModel_Delete_AuxFMU (void *MP)
{
    /* Do not change */
    tAuxData *mp = (struct tAuxData*)MP;

    /* Free strorage for Head */
    free(&mp->Head);

    /* Replace this block with user model if FMU other than MySuspKnC_FMU is used */
    {
	/* Free storage for user model */
	MyModel_Delete(mp->MyAuxModel);
    }
}

int
Susp_KnC_Register_MyModel_AuxFMU (void)
{
    tModelClassDescr m;

    memset (&m, 0, sizeof(m));
    m.SuspKnC.VersionId =	ThisVersionId;
    m.SuspKnC.New =		MyModel_New_AuxFMU;
    m.SuspKnC.Calc =		MyModel_Calc_AuxFMU;
    m.SuspKnC.Delete =		MyModel_Delete_AuxFMU;
    /* Should only be used if the model doesn't read params from extra files */
    m.SuspKnC.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_SuspKnC, ThisModelKindAuxFMU, &m);
}
