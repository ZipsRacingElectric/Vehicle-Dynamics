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
** Simple motor Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    Motor_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "Powertrain.Motor";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


struct tMyModel {
    double Inert;
    double rotv;
    double Ratio;
    double Trq;
};

static void MyModel_Delete (void *MP);


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
MyModel_New (struct tInfos *Inf, struct tPTMotorCfgIF *CfgIF,
	     const char *KindKey , const char *Ident)
{
    struct tMyModel *mp = NULL;
    const char *ModelKind;
    char MsgPre[64];
    int VersionId = 0;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_PTMotor, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*)calloc(1,sizeof(*mp));

    /* get CfgIF parameters */
    if (Motor_GetCfgOutIF (Inf, CfgIF, ModelKind, Ident) != 0)
	goto ErrorReturn;

    /* CfgIF output -> Model */
    mp->Ratio  = CfgIF->Ratio;

    mp->Inert  = fabs(iGetDbl (Inf, "PowerTrain.MotorISG.I"));
    mp->Inert *= mp->Ratio * mp->Ratio;

    /* CfgIF output: verification if the parametrization corresponds to the model */
    if (CfgIF->TrqMot_max==NULL) {
	LogErrF (EC_Init, "%s: no characteristic for motor torque", MsgPre);
	goto ErrorReturn;
    }
    //CfgIF->TrqGen_max in the model ignored

    return mp;

    ErrorReturn:
	MyModel_Delete (mp);
	return NULL;
}

static int
MyModel_Calc (void *MP, struct tPTMotorIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    const tPTMotorCfgIF *CfgIF = IF->CfgIF;

    mp->rotv = IF->rotv * mp->Ratio;
    mp->Trq  = LMEval (CfgIF->TrqMot_max, fabs(mp->rotv)) * IF->Load;

    IF->Trq     = mp->Trq * mp->Ratio;
    IF->PwrElec = IF->Trq * IF->rotv;
    IF->Inert   = mp->Inert;

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
Motor_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.PTMotor.VersionId =	ThisVersionId;
    m.PTMotor.New =		MyModel_New;
    m.PTMotor.Calc =		MyModel_Calc;
    m.PTMotor.Delete =		MyModel_Delete;
    m.PTMotor.DeclQuants =	MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.PTMotor.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_PTMotor, ThisModelKind, &m);
}
