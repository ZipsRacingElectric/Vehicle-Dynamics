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
** Simple clutch Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    Clutch_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "PowerTrain.Clutch";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


struct tMyModel {
    /* Parameters */
    double	c;	/* Torque coefficient [Nms/rad] */
    double	I_in;	/* Inertia in and out [kgm^2] */
    double	I_out;
    double	Trq_max;

    /* Variables */
    double	drotv;	/* Difference of rotation velocity [rad/s] */
    double	Trq;	/* Clutch torque [Nm] */
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
MyModel_New (tInfos *Inf, tPTClutchCfgIF *CfgIF, const char *KindKey, const char *Ident)
{
    struct tMyModel *mp = NULL;
    const char *ModelKind, *key;
    char MsgPre[64];
    int VersionId = 0;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_PTClutch, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*)calloc(1,sizeof(*mp));

    /* get CfgIF parameters */
    if (Clutch_GetCfgOutIF (Inf, CfgIF, ModelKind, Ident) != 0)
	goto ErrorReturn;

    /* Inertia In/Out */
    mp->I_in    = iGetDblOpt(Inf, "PowerTrain.Clutch.I_in",     1e-4);
    mp->I_out   = iGetDblOpt(Inf, "PowerTrain.Clutch.I_out",    1e-4);

    mp->Trq_max = iGetDblOpt(Inf, "PowerTrain.Clutch.Trq_max",  500.0);

    key = "PowerTrain.MyClutch.c";
    mp->c = iGetDbl(Inf, key);
    if (mp->c <= 0.0) {
	LogErrF(EC_Init,
	    "%s: torque coefficient '%s' must be positive and non zero", MsgPre, key);
	goto ErrorReturn;
    }

    /* CfgIF output: verification if the parametrization corresponds to the model */
    if (CfgIF->ClKind != ClKind_Friction) {
	LogErrF (EC_Init, "%s: model supports only a friction clutch", MsgPre);
	goto ErrorReturn;
    }

    return mp;

    ErrorReturn:
	free (mp);
	return NULL;
}


static int
MyModel_PreSimSetup (void *MP, struct tPTClutchPreSimIF *PreSimIF,
		      	       struct tPTClutchIF       *IF)
{
    IF->rotv_in = PreSimIF->rotv_in;

    return 0;
}


static int
MyModel_Calc (void *MP, struct tPTClutchIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    double rota_in;

    /* Clutch Friction Torque */
    mp->drotv = IF->rotv_in - IF->rotv_out;
    mp->Trq   = mp->drotv * mp->c;
    M_BOUND_ABS (mp->Trq_max, mp->Trq);
    mp->Trq  *= M_BOUND (0.0, 1.0, 1.0 - IF->Pos);

    /* Integration of Clutch DOF */
    rota_in      = (IF->Trq_in - mp->Trq) / (IF->Inert_in + mp->I_in);
    IF->rotv_in += rota_in     * dt;
    IF->rot_in  += IF->rotv_in * dt;

    IF->Trq_SuppInert = rota_in * (IF->Inert_in + mp->I_in);

    IF->Trq_out   = mp->Trq;
    IF->Inert_out = mp->I_out;

    IF->i_TrqIn2Out = 1.0 - IF->Pos;

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


static int
MyModel_ModelCheck (void *MP, struct tInfos *Inf)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    FILE *fp;

    if ((fp = ModelCheck_GetDesignFile()) == NULL)
	return 0;

    fprintf (fp, "### Clutch.Kind = %s\n", ThisModelKind);
    fprintf (fp, "Clutch.I_in =             %10.7f\n", mp->I_in);
    fprintf (fp, "Clutch.I_out =            %10.7f\n", mp->I_out);
    fprintf (fp, "\n");

    return 0;
}


int 
Clutch_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.PTClutch.VersionId =	ThisVersionId;
    m.PTClutch.New =		MyModel_New;
    m.PTClutch.Calc =		MyModel_Calc;
    m.PTClutch.Delete =		MyModel_Delete;
    m.PTClutch.DeclQuants =	MyModel_DeclQuants;
    m.PTClutch.ModelCheck =	MyModel_ModelCheck;
    m.PTClutch.PreSimSetup =	MyModel_PreSimSetup;
    /* Should only be used if the model doesn't read params from extra files */
    m.PTClutch.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_PTClutch, ThisModelKind, &m);
}
