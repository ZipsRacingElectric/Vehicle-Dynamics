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
** Simple aerodynamics Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    Aero_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "Aero";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


struct tMyModel {
    /* Parameters */
    double	Area;		/* Reference area [m^2] */
    double	Length;		/* Reference length [m] */
    double	PoA_1[3];	/* Point of attack [m] */
    tLM		*cMap;		/* cX coefficients map */
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



static void *
MyModel_New (tInfos *Inf, const char *KindKey)
{
    struct tMyModel *mp = NULL;
    int n;
    const char *key;
    char MsgPre[64];
    double dvec[512];

    const char *ModelKind;
    int VersionId = 0;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_Aero, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*)calloc(1,sizeof(*mp));

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    mp->Area =   iGetDbl (Inf, "MyAero.Area");
    mp->Length = iGetDbl (Inf, "MyAero.Length");
    iGetTable (Inf, "MyAero.PoA_1", mp->PoA_1, 3, 1, &n);

    key = "MyAero.Coeff";
    if ((iGetTable(Inf, key, dvec, 512, 7, &n)) < 0) {
	LogErrF (EC_Init, "%s: Error while reading '%s'", MsgPre, key);
	goto ErrorReturn;
    } else {
	int i, k;
	tMData *md = MData3DNew(n, 1, 0, 6);
	
	for (i=0; i < n; i++) {
	    MDataX1SetVal (md, i, dvec[i]*deg2rad);
	    for (k=0; k < 6; k++) {
		MData3DSetVal (md, i, 0, 0, k, dvec[i+(k+1)*n]);
	    }
	}
	mp->cMap = LMInitMD(md, ArrayNotAequi);
	MDataDelete (md);
	md = NULL;
	
	if (mp->cMap == NULL) {
	    LogErrF(EC_Init, "%s: Can't generate aero mapping", MsgPre);
	    goto ErrorReturn;
	}
    }

    return mp;

  ErrorReturn:
    if (mp->cMap != NULL)
	LMDelete (mp->cMap);
    free (mp);
    return NULL;
}


static int
MyModel_Calc (void *MP, tAeroIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    double vel, x, cX[6];

    vel = sqrt(IF->ApproachVel_1[0]*IF->ApproachVel_1[0]
	       + IF->ApproachVel_1[1]*IF->ApproachVel_1[1]);

    /* PoA_1 */
    IF->PoA_1[0] = mp->PoA_1[0];
    IF->PoA_1[1] = mp->PoA_1[1];
    IF->PoA_1[2] = mp->PoA_1[2];

    /* cX coefficients depending on tau_1 [rad] */
    LMVEval (mp->cMap, IF->tau_1, cX);

    /* Fr1_1 */
    x = Env.AirDensity * 0.5 * mp->Area * vel*vel;
    IF->Frc_1[0] = x * cX[0];
    IF->Frc_1[1] = x * cX[1];
    IF->Frc_1[2] = x * cX[2];

    /* Trq_1 */
    x = Env.AirDensity * 0.5 * mp->Area * mp->Length * vel*vel;
    IF->Trq_1[0] = x * cX[3];
    IF->Trq_1[1] = x * cX[4];
    IF->Trq_1[2] = x * cX[5];

    return 0;
}


static void
MyModel_Delete (void *MP)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Park the dict entries for dynamically allocated variables before deleting */
    MyModel_DeclQuants_dyn (mp, 1);

    if (mp->cMap != NULL)
	LMDelete (mp->cMap);

    free(mp);
}


int 
Aero_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.Aero.VersionId =		ThisVersionId;
    m.Aero.New =		MyModel_New;
    m.Aero.Calc =		MyModel_Calc;
    m.Aero.Delete =		MyModel_Delete;
    m.Aero.DeclQuants =		MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.Aero.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_Aero, ThisModelKind, &m);
}
