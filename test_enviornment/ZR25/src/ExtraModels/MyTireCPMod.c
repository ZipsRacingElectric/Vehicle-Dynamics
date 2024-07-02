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
** Simple tire Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    TireCPMod_Register_MyModel ();
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
    /* Model specific Parameters */
    int TireNo;
    double HitTime;
    double HitDuration;
    double z_offset;
};


static void
MyModel_DeclQuants_dyn (struct tMyModel *mp, int park)
{
    tDDefault *df;
    static struct tMyModel MyModel_Dummy = {0};
    if (park)
	mp = &MyModel_Dummy;

    /* Define here dict entries for dynamically allocated variables. */

    df = DDefaultCreate ("MyTireCPMod.");
    DDefDouble4 (df, "HitTime",		"s",	&mp->HitTime,  DVA_VC);
    DDefDouble4 (df, "z_offset",	"m",	&mp->z_offset, DVA_VC);
    DDefaultDelete (df);
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
MyModel_New (
    struct tInfos	*inf,
    int			TireNo,
    int			Side,
    int			Twin)
{
    struct tMyModel *mp = (struct tMyModel*)calloc(1,sizeof(*mp));
    mp->TireNo = TireNo;

    mp->HitTime      = iGetDblOpt (inf, "MyTireCPMod.HitTime",     10.0);
    mp->HitDuration  = iGetDblOpt (inf, "MyTireCPMod.HitDuration", 0.01);
    mp->z_offset     = iGetDblOpt (inf, "MyTireCPMod.z_offset",    0.2);

    return mp;
}


static int
MyModel_Calc (void *MP, tTireCPModIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Hit on the both rear tires at time=TimeHit by an additional contact point z offset */
    if ((mp->TireNo==2 || mp->TireNo==3) && SimCore.Time>mp->HitTime
	 && SimCore.Time<mp->HitTime+mp->HitDuration) {
	IF->CP_0[2] += mp->z_offset;
	Log ("T=%g hit on tire=%d\n", SimCore.Time, mp->TireNo);
    }

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
TireCPMod_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.TireCPMod.VersionId =	ThisVersionId;
    m.TireCPMod.New =		MyModel_New;
    m.TireCPMod.Calc =		MyModel_Calc;
    m.TireCPMod.Delete =	MyModel_Delete;
    m.TireCPMod.DeclQuants =	MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.TireCPMod.ParamsChanged = ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_TireCPMod, ThisModelKind, &m);
}


/* TireCPMod for Four Poster **************************************************/

static const char ThisModelKind_FP[]  = "MyModel_FourPoster";
static const int  ThisVersionId_FP    = 1;


struct tMyModel_FP {
    /* Model specific Parameters */
    int TireNo;
    int Twin;
    double z_offset;
};


static void
MyModel_DeclQuants_dyn_FP (struct tMyModel_FP *mp, int park)
{
    tDDefault *df;
    static struct tMyModel_FP MyModel_Dummy_FP = {0};

    if (park) {
	MyModel_Dummy_FP.Twin = mp->Twin;
	MyModel_Dummy_FP.TireNo = mp->TireNo;
	mp = &MyModel_Dummy_FP;	
    }

    const char *buft = mp->Twin ? ".Twin" : "";

    /* Define here dict entries for dynamically allocated variables. */
    df = DDefaultCreate(NULL);

    if (mp->TireNo >= TireNo_FirstTrailer2Tire)
	DDefPrefix(df, "MyTireCPMod_FP.Tr2.%s%s.", Tire_TireNo_Str(mp->TireNo), buft);
    else if (mp->TireNo >= TireNo_FirstTrailerTire)
	DDefPrefix(df, "MyTireCPMod_FP.Tr.%s%s.", Tire_TireNo_Str(mp->TireNo), buft);
    else
	DDefPrefix(df, "MyTireCPMod_FP.%s%s.", Tire_TireNo_Str(mp->TireNo), buft);

    DDefDouble4 (df, "z_offset",	"m",	&mp->z_offset, DVA_VC);
    DDefaultDelete (df);
}


static void
MyModel_DeclQuants_FP (void *MP)
{
    struct tMyModel_FP *mp = (struct tMyModel_FP *)MP;

    if (mp == NULL) {
	/* Define here dict entries for non-dynamically allocated (static) variables. */

    } else {
	MyModel_DeclQuants_dyn_FP (mp, 0);
    }
}


static void *
MyModel_New_FP (struct tInfos *inf, int TireNo, int Side, int Twin)
{
    struct tMyModel_FP *mp = (struct tMyModel_FP*)calloc(1,sizeof(*mp));
    mp->TireNo = TireNo;

    return mp;
}


static int
MyModel_Calc_FP (void *MP, tTireCPModIF *IF, double dt)
{
    struct tMyModel_FP *mp = (struct tMyModel_FP *)MP;

    IF->CP_0[2] += mp->z_offset;

    return 0;
}


static void
MyModel_Delete_FP (void *MP)
{
    struct tMyModel_FP *mp = (struct tMyModel_FP *)MP;

    /* Park the dict entries for dynamically allocated variables before deleting */
    MyModel_DeclQuants_dyn_FP (mp, 1);
    free (mp);
}


int 
TireCPMod_Register_MyModelFourPoster (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.TireCPMod.VersionId =	ThisVersionId_FP;
    m.TireCPMod.New =		MyModel_New_FP;
    m.TireCPMod.Calc =		MyModel_Calc_FP;
    m.TireCPMod.Delete =	MyModel_Delete_FP;
    m.TireCPMod.DeclQuants =	MyModel_DeclQuants_FP;
    /* Should only be used if the model doesn't read params from extra files */
    m.TireCPMod.ParamsChanged = ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_TireCPMod, ThisModelKind_FP, &m);
}
