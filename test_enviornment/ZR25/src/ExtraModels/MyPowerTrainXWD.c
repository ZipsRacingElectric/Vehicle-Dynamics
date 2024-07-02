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
**    PowerTrainXWD_Register_MyModel ();
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

static const char ThisModelClass[] = "PowerTrain.PTXWD";
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
MyModel_New (tInfos *Inf, struct tPowerTrainXWD_CfgIF *CfgIF, const char *KindKey)
{
    struct tMyModel *mp = NULL;
    int VersionId = 0;
    char MsgPre[64];
    const char *ModelKind;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_PowerTrainXWD, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*)calloc(1,sizeof(*mp));

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    /* get CfgIF parameters */
    if (PowerTrainXWD_GetCfgOutIF (Inf, CfgIF, ModelKind) != 0)
	goto ErrorReturn;

    if (CfgIF->nWheels != NWHEEL) {
	LogErrF(EC_Init, "%s: model supports only a four wheel vehicle", MsgPre);
	goto ErrorReturn;
    }

    /* CfgIF output: verification if the parametrization corresponds to the model */
    if (CfgIF->DriveLine.iDiff_mean <= 0) {
	LogErrF (EC_Init, "%s: mean driveline ratio must be positive and non zero", MsgPre);
	goto ErrorReturn;
    }
    if (CfgIF->DriveLine.DriveSourcePos[0] != Diff_Front ||
	CfgIF->DriveLine.DriveSourcePos[1] != NoPosition  ||
        CfgIF->DriveLine.DriveSourcePos[2] != NoPosition  ||
        CfgIF->DriveLine.DriveSourcePos[3] != NoPosition) {
	LogErrF (EC_Init, "%s: model supports only one drive source at center differential", MsgPre);
	goto ErrorReturn;
    }
    if (CfgIF->PTKind != PTKind_BEV) {
	LogErrF (EC_Init, "%s: model supports only electrical powertrain", MsgPre);
	goto ErrorReturn;
    }
    if (CfgIF->GearBox.GBKind != GBKind_NoGearBox) {
	LogErrF (EC_Init, "%s: model supports only no gearbox", MsgPre);
	goto ErrorReturn;
    }
    if (CfgIF->nMotor != 1) {
	LogErrF (EC_Init, "%s: model supports only one electrical motor", MsgPre);
	goto ErrorReturn;
    }

    return mp;

    ErrorReturn:
	free (mp);
	return NULL;
}


static int
MyModel_Calc (void *MP, struct tPowerTrainXWD_IF *IF, double dt)
{
    const tPowerTrainXWD_CfgIF *CfgIF            =  IF->CfgIF;
#ifdef __cplusplus
    struct tPowerTrainXWD_IF::tPowerTrainXWD_IF_MotorIn  *MIn    = &IF->MotorIn[0];
    struct tPowerTrainXWD_IF::tPowerTrainXWD_IF_MotorOut *MOut   = &IF->MotorOut[0];
    struct tPowerTrainXWD_IF::tPowerTrainXWD_IF_BattOut  *BattHV = &IF->BattHVOut;
#else
    struct tPowerTrainXWD_IF_MotorIn  *MIn       = &IF->MotorIn[0];
    struct tPowerTrainXWD_IF_MotorOut *MOut      = &IF->MotorOut[0];
    struct tPowerTrainXWD_IF_BattOut  *BattHV    = &IF->BattHVOut;
#endif
    int i;

    /** PowerSupply **/
    /* in this example, provide unlimited Power to the electrical motor */
    BattHV->Pwr_max = 400;

    /** electrical Motor **/
    MOut->rotv = (IF->WheelIn[0].rotv + IF->WheelIn[1].rotv) / 2 * CfgIF->DriveLine.iDiff_mean;
    MOut->Trq = LMEval (CfgIF->Motor[0].TrqMot_max, fabs(MOut->rotv)) * MIn->Load;
    MOut->PwrElec = MOut->Trq * MOut->rotv;

    /** DriveLine **/
    /* Calculate the drive torque at wheels */
    for (i=0; i<2; i++)
	IF->WheelOut[i].Trq_Drive = MOut->Trq * 0.5 * CfgIF->DriveLine.iDiff_mean;

    /* Support torque to the wheel carrier */
    for (i=0; i<NWHEEL; i++)
	IF->WheelOut[i].Trq_Supp2WC = -IF->WheelOut[i].Trq_Drive;

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
PowerTrainXWD_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.PowerTrainXWD.VersionId =		ThisVersionId;
    m.PowerTrainXWD.New =		MyModel_New;
    m.PowerTrainXWD.Calc =		MyModel_Calc;
    m.PowerTrainXWD.Delete =		MyModel_Delete;
    m.PowerTrainXWD.DeclQuants =	MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.PowerTrainXWD.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_PowerTrainXWD, ThisModelKind, &m);
}
