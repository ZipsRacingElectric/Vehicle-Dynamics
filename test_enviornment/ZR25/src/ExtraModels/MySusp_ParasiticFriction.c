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
**  suspension Model parasitic friction: LuGre model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    int Susp_ParasiticFriction_Register_MyModel (void);
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Car.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"
#include "Log.h"

static const char ThisModelClass[] = "Susp_ParasiticEffects";
static const char ThisModelKind[]  = "MyFrictionModel";
static const int  ThisVersionId    = 1;

/* model parameters */
struct tMyModel {
    int SuspModID;
    struct tMyFrcParasiticFriction {
	double z_dot;
	double z;
	double sigma0;
	double sigma1;
	double v0;
	double mu_c;
	double PreloadForce;
    } ParasiticFriction;
};


static void
MyModel_DeclQuants_dyn (struct tMyModel *mp, int park, int SuspModID)
{
    static struct tMyModel MyModel_Dummy;
    memset (&MyModel_Dummy, 0, sizeof(struct tMyModel));
    if (park)
	mp = &MyModel_Dummy;

    char Quant[40];
    char MdlPre[4];

    if (SuspModID >= SuspMod_FirstTrailer2SId)
	sprintf(MdlPre, "Tr2");
    else if (SuspModID >= SuspMod_FirstTrailerSId)
	sprintf(MdlPre, "Tr");
    else
	sprintf(MdlPre, "Car");

    const char *s = SuspMod_StrId_uc(SuspModID);

    sprintf(Quant, "%s.Parasitic.Friction%s.z", MdlPre, s);
    DDefDouble(NULL, Quant, "m", &mp->ParasiticFriction.z, DVA_None);
    sprintf(Quant, "%s.Parasitic.Friction%s.z_dot", MdlPre, s);
    DDefDouble(NULL, Quant, "m/s", &mp->ParasiticFriction.z_dot, DVA_None);
}


static void
MyModel_DeclQuants (void *MP, int SuspModID)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    if (mp == NULL) {
	/* Define here dict entries for non-dynamically allocated (static) variables. */

    } else {
	MyModel_DeclQuants_dyn (mp, 0, SuspModID);
    }
}


static void *
MyModel_New (struct tInfos *Inf, struct tSusp_ParasiticEffectsCfgIF *CfgIF, const char *KindKey)
{
    struct tMyModel *mp = NULL;
    int VersionId = 0;
    char MsgPre[64];
    const char *ModelKind;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_Susp_ParasiticEffects, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    mp = (struct tMyModel*)calloc(1, sizeof(*mp));

    char Key[50];
    const char *s = SuspMod_SuspPre(CfgIF->SuspModID);

    mp->SuspModID = CfgIF->SuspModID;
    
    /* LuGre */
    sprintf(Key, "SFH.%s.Parasitic.Friction.LuGre.Stiffness", s);
    mp->ParasiticFriction.sigma0 = iGetDbl(Inf, Key);
    sprintf(Key, "SFH.%s.Parasitic.Friction.LuGre.Damping", s);
    mp->ParasiticFriction.sigma1 = iGetDbl(Inf, Key);
    sprintf(Key,"SFH.%s.Parasitic.Friction.LuGre.v0", s);
    mp->ParasiticFriction.v0 = iGetDbl(Inf, Key);
    sprintf(Key, "SFH.%s.Parasitic.Friction.LuGre.mu_c", s);
    mp->ParasiticFriction.mu_c = iGetDbl(Inf, Key);
    sprintf(Key, "SFH.%s.Parasitic.Friction.LuGre.PreloadForce", s);
    mp->ParasiticFriction.PreloadForce = iGetDbl(Inf, Key);

    return mp;
}

//Dir function
static double
dir (double v, double veps)
{
    double veps_rezip = 1.0 / veps;

    if (fabs(v) >= veps){
	return v/fabs(v);
    } else {
	return 	v * veps_rezip *(1.5 * fabs(v) * veps_rezip
				 - 0.5 * (fabs(v) * veps_rezip) * (fabs(v) * veps_rezip) * (fabs(v) * veps_rezip));
    }
}

static int
MyModel_Calc (void *MP, tSusp_ParasiticEffectsIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel*)MP;

    double g;

    g = exp( - (IF->vWheelCarrier_z * IF->vWheelCarrier_z)/(mp->ParasiticFriction.v0 * mp->ParasiticFriction.v0));
    mp->ParasiticFriction.z_dot =
    	g * IF->vWheelCarrier_z
    	+ (1. - g) * (    1. / mp->ParasiticFriction.sigma1 * mp->ParasiticFriction.mu_c * dir(IF->vWheelCarrier_z, 0.01)
    			- mp->ParasiticFriction.sigma0 / mp->ParasiticFriction.sigma1 * mp->ParasiticFriction.z);
    mp->ParasiticFriction.z += mp->ParasiticFriction.z_dot * dt;

    IF->Force = mp->ParasiticFriction.PreloadForce *   (  mp->ParasiticFriction.sigma0 * mp->ParasiticFriction.z
    							+ mp->ParasiticFriction.sigma1 * mp->ParasiticFriction.z_dot);

    return 0;
}


static void
MyModel_Delete (void *MP)
{
    struct tMyModel *mp = (struct tMyModel*)MP;

    if (mp == NULL)
        return;

    /* Park the dict entries for dynamically allocated variables before deleting */
    MyModel_DeclQuants_dyn (mp, 1, mp->SuspModID);

    free(mp);
    mp = NULL;
}


int
Susp_ParasiticFriction_Register_MyModel (void)
{
    tModelClassDescr m;

    memset (&m, 0, sizeof(m));
    m.Susp_ParasiticEffects.VersionId =	    	ThisVersionId;
    m.Susp_ParasiticEffects.New =		MyModel_New;
    m.Susp_ParasiticEffects.Calc =		MyModel_Calc;
    m.Susp_ParasiticEffects.Delete =	    	MyModel_Delete;
    m.Susp_ParasiticEffects.DeclQuants =	MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.Susp_ParasiticEffects.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_Susp_ParasiticEffects, ThisModelKind, &m);
}
