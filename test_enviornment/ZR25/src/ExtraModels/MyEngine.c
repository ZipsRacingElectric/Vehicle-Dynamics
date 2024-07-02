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
** Simple engine Model with two look-up tables
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    Engine_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "PowerTrain.Engine";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


struct tMyModel {

    double	TrqKl15Off;	/* Torque for Kl15=0                     [Nm] */
    double	rotvIdle;	/* Engine idle speed		      [rad/s] */
    double	rotvOff;	/* Engine off speed		      [rad/s] */
    double	rotvMin;	/* Engine minimum speed		      [rad/s] */
    double	rotvMax;	/* Engine maximum speed		      [rad/s] */
    tLM		*TrqFull;	/* Look-Up Table 1D for full-load torque [Nm] */
    tLM		*TrqDrag;	/* Look-Up Table 1D for drag torque      [Nm] */
    double	I_out;		/* Inertia outshaft                   [kgm^2] */
};

static void MyModel_Delete (void *MP);


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


/* Model output parameters in the configuration struct CfgIF, which are required
   by CarMaker, are read in before the MyModel_New() function.
   - The parametrization of these parameters is supported by the GUI.
   - These output parameters can be used internally by the model in same way like
     the input parameters
*/
static void *
MyModel_New (
    struct tInfos	  *Inf,
    struct tPTEngineCfgIF *CfgIF,
    const char		  *KindKey,
    const char		  *Ident)
{
    struct tMyModel *mp = NULL;
    char 	MsgPre[64];
    int		VersionId = 0;
    const char *ModelKind;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_PTEngine, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*)calloc(1,sizeof(*mp));

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    /* get CfgIF parameters */
    if (Engine_GetCfgOutIF (Inf, CfgIF, ModelKind) != 0)
	goto ErrorReturn;

    /* Get parameters from file */
    mp->TrqKl15Off = iGetDbl(Inf, "PowerTrain.MyEngine.TrqKl15Off");

    /* Inertia */
    mp->I_out = iGetDblOpt (Inf, "PowerTrain.Engine.I", 0.1);

    /* CfgIF output -> Model */
    mp->TrqFull = CfgIF->TrqFull;
    mp->TrqDrag = CfgIF->TrqDrag;

    /* CfgIF output: verification if the parametrization corresponds to the model */
    if (mp->TrqFull==NULL || mp->TrqDrag==NULL) {
	LogErrF (EC_Init, "%s: missing engine torque characteristic", MsgPre);
	goto ErrorReturn;
    }

    return mp;

  ErrorReturn:
    MyModel_Delete (mp);
    return NULL;
}


/* Torque zero for rotation under reference speed */
static double
Fac4VelZero (double vel, double vel_ref)
{
    double absvel = fabs(vel);

    if (absvel >= vel_ref)
	return 1.0;
    else
	return 0.5 * (1.0 - cos(M_PI * absvel / vel_ref));
}


static int
MyModel_Calc (void *MP, struct tPTEngineIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    double TrqFull, TrqDrag, Load = M_BOUND(0.0, 1.0, IF->Load);

    if (IF->FuelLevel==0.0)
	Load = 0.0;

    if (IF->Ignition) {
	double rotv = M_MAX(IF->rotv, mp->rotvMin);
	
	/* Use the torque from the Look-Up tables = f(rotv, Load) */
	TrqFull = LMEval(mp->TrqFull, rotv);
	TrqDrag = LMEval(mp->TrqDrag, rotv);

	IF->Trq = TrqDrag + Load * (TrqFull-TrqDrag);
    } else {
	/* Use the torque for Kl15=0 */
	IF->Trq = mp->TrqKl15Off;
    }
    IF->Trq  *= M_SGN(IF->rotv) * Fac4VelZero(IF->rotv, 2.0*rpm2radsec);

    IF->Inert = mp->I_out;

    return 0;
}


static int
MyModel_ModelCheck (void *MP, struct tInfos *Inf)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    FILE *fp;

    if ((fp = ModelCheck_GetDesignFile()) == NULL)
	return 0;

    fprintf (fp, "### Engine.Kind = %s\n", ThisModelKind);
    fprintf (fp, "Engine.I =                %10.7f\n", mp->I_out);
    fprintf (fp, "\n");

    return 0;
}


static void
MyModel_Delete (void *MP)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Park the dict entries for dynamically allocated variables before deleting */
    MyModel_DeclQuants_dyn (mp, 1);

    free(mp);
}


int 
Engine_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.PTEngine.VersionId =	ThisVersionId;
    m.PTEngine.New =		MyModel_New;
    m.PTEngine.Calc =		MyModel_Calc;
    m.PTEngine.Delete =		MyModel_Delete;
    m.PTEngine.DeclQuants =	MyModel_DeclQuants;
    m.PTEngine.ModelCheck =	MyModel_ModelCheck;
    /* Should only be used if the model doesn't read params from extra files */
    m.PTEngine.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_PTEngine, ThisModelKind, &m);
}
