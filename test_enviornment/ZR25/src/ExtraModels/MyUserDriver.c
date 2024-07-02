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
** Simple UserDriver for longitudinal dynamics to reach and maintain a given
** target velocity. It is not capable of starting or changing gears.
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    UserDriver_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "Vehicle.h"
#include "MyModels.h"

static const char ThisModelClass[] = "UserDriver";
static const char ThisModelKind[]  = "MyUserDriver";
static const int  ThisVersionId    = 1;


struct tMyModel {
    double	vTrgt;		/* Target velocity */
    double 	tPrev;		/* Preview time */
    double	kp, ki, kd;	/* PID-Controller */
    double	CtrlInt;	/* Integrator */
    double	CtrlDiff;	/* Differentiator */
    double 	TLPF, f_old;	/* LPF */
    int 	Init;		/* Initialisation */
};


static void
MyModel_DeclQuants_dyn (struct tMyModel *mp, int park)
{
    static struct tMyModel MyModel_Dummy = {0};
    if (park)
	mp = &MyModel_Dummy;

    /* Define here dict entries for dynamically allocated variables. */
}

/*
** MyModel_DeclQuants ()
**
** Defining DataDictionary Quantities
**
** Call:
** - at the beginning of every TestRun, after the call for the _New function
*/

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

/*
** MyModel_New ()
**
** Initialising the model
**
** Call:
** - one times at the beginning of every TestRun
*/

static void *
MyModel_New (tInfos *Inf, const char *KindKey)
{
    struct tMyModel *mp = NULL;
    const char *ModelKind;
    char MsgPre[64];
    int VersionId = 0;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_UserDriver, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    if ((mp = (struct tMyModel*)calloc(1,sizeof(*mp))) == NULL)
        return NULL;

    /* Parameters */
    mp->kp = iGetDblOpt(Inf, "UserDriver.kp", 0.5);
    mp->ki = iGetDblOpt(Inf, "UserDriver.ki", 0.5);
    mp->kd = iGetDblOpt(Inf, "UserDriver.kd", 0.0);
    mp->tPrev = iGetDblOpt(Inf, "UserDriver.tPreview", 5.0);
    mp->vTrgt = iGetDblOpt(Inf, "UserDriver.TargetVel", 70.0) / 3.6;
    mp->TLPF = iGetDblOpt(Inf, "UserDriver.TLPF", 0.5);
    mp->f_old = mp->CtrlInt = mp->CtrlDiff =  0.0;
    mp->Init = 1;

    return mp;
}

static double
Filter_PT1(double delta, double dt, struct tMyModel *mp)
{
    return mp->f_old + (delta - mp->f_old) * dt/(mp->TLPF + dt);
}


static double
PIDCtrl_Calc(double delta, double dt, struct tMyModel *mp) {
    double ret;

    if (mp->Init) {
        mp->CtrlDiff = delta;
	mp->Init = 0;
    }

    /* Low-pass filter */
    delta = mp->f_old = Filter_PT1(delta, dt, mp);

    ret = mp->kp * delta;			 /* P */
    ret += mp->ki * mp->CtrlInt * dt;		 /* I */
    ret += mp->kd * (delta - mp->CtrlDiff) / dt; /* D */

    /* Limitation */
    if (ret > 1)
	ret = 1;
    else if (ret < -1)
	ret = -1;
    else
	mp->CtrlInt += delta; /* Anti-Windup ("Clamping") */

    mp->CtrlDiff = delta;

    return ret;
}


/*
** MyModel_Calc ()
**
** Calculation of the model
**
** Call:
** - every cycle
*/

static int
MyModel_Calc (void *MP, struct tUserDriverIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    double var, delta, vProg;

    vProg = Vehicle.PoI_Vel_1[0] + Vehicle.PoI_Acc_1[0] * mp->tPrev;
    delta = mp->vTrgt - vProg;

    /* PID-Controller */
    var = PIDCtrl_Calc(delta, dt, mp);

    /* Gas or Brake */
    if (var >= 0){
	IF->Out.Gas   = IF->In.Clutch >= 0.1 ? 0 : var;
	IF->Out.Brake = 0;
    } else {
	IF->Out.Gas   = 0;
	IF->Out.Brake = -var;
    }

    return 0;
}

/*
** MyModel_Delete ()
**
** Uninitialising the model
**
** Call:
** - at the end of every TestRun
*/

static void
MyModel_Delete (void *MP)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Park the dict entries for dynamically allocated variables before deleting */
    MyModel_DeclQuants_dyn (mp, 1);
    free (mp);
}


int 
UserDriver_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.UserDriver.VersionId =	ThisVersionId;
    m.UserDriver.New =		MyModel_New;
    m.UserDriver.Calc =		MyModel_Calc;
    m.UserDriver.Delete =	MyModel_Delete;
    m.UserDriver.DeclQuants =	MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    // m.UserDriver.ParamsChanged = ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_UserDriver, ThisModelKind, &m);
}
