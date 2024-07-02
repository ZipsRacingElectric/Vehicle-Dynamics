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
*/
#include <Global.h>

#if defined(WIN32)
#  include <windows.h>
#endif
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "Log.h"
#include "DrivMan.h"
#include "DataDict.h"
#include "InfoUtils.h"
#include "SimCore.h"
#include "ModelManager.h"
#include "MathUtils.h"

#include "Vehicle/Sensor_Object.h"
#include "Vehicle.h"
#include "VehicleControl.h"
#include "Vehicle/VehicleControlApps.h"

#define NOTSET	-99999

tAccelCtrl AccelCtrl;

static double c_i;
static tDesrAccelFunc User_DesrAccelFunc = NULL;


static void
AccelCtrl_DeclQuants (void *MP)
{
    tACC_ECU *ACC_ECU = &AccelCtrl.ACC_ECU;
    tDDictEntry *e;
    tDDefault   *df;

    if (MP==NULL)
	return;

    df = DDefaultCreate("AccelCtrl.");

    DDefDouble4 (df, "DesiredAx",          "m/s^2", &AccelCtrl.DesrAx,        DVA_DM);

    e = DDefChar (df, "ACC.IsActive",      "",      &ACC_ECU->IsActive,       DVA_DM);
    DDefStates(e, 2, 0);

    DDefDouble4 (df, "ACC.DesiredSpd",     "m/s",   &ACC_ECU->DesrSpd,        DVA_DM);
    DDefDouble4 (df, "ACC.DesiredDist",	   "m",     &ACC_ECU->DesrDist,       DVA_None);
    DDefDouble4 (df, "ACC.DesiredTGap",	   "s",     &ACC_ECU->DesrTGap,       DVA_DM);
    DDefDouble4 (df, "ACC.DesiredAx",	   "m/s^2", &ACC_ECU->DesrAx,         DVA_DM);
    DDefDouble4 (df, "ACC.Time2Collision", "s",     &ACC_ECU->Time2Collision, DVA_None);

    DDefaultDelete(df);
}


static void
DesrAccelFunc_ACC (double dt)
{
    double ax, ax_sc, delta_ds;
    tACC_ECU *ACC_ECU = &AccelCtrl.ACC_ECU;
    struct tObjectSensor *s = &ObjectSensor[ACC_ECU->RefObjectSensorId];

    /* Driver Brake Limit */
    if (DrivMan.Brake > ACC_ECU->BrakeThreshold)
	ACC_ECU->IsActive = 0;

    /* Time until collision */
    if (s->relvTarget.NearPnt.dv_p < 0) {
	ACC_ECU->Time2Collision =  s->relvTarget.NearPnt.ds_p
			       / -s->relvTarget.NearPnt.dv_p;
    } else {
	ACC_ECU->Time2Collision = 0;
    }

    if (!ACC_ECU->IsActive) {
	/* ACC off -> set desired speed to current car speed */
	ACC_ECU->DesrSpd = Vehicle.v;
	ACC_ECU->DesrAx  = NOTSET;
	goto SetAccelCtrl;
    }

    /* ACC active */
    if (s->Targ_Dtct) {
	/* if target detected set desired distance,
       DesrDistance[m] = Target.v[m/s] * 3.6 / Desired Time Gap(Init= 1.8[s])
       or if target stand still DSMIN: 20[m] distance */
	ACC_ECU->DesrDist = M_MAX( ((Vehicle.v +
		s->relvTarget.NearPnt.dv_p) * ACC_ECU->DesrTGap), ACC_ECU->dsmin);
	
	/* Distance Control Algorithm: result = desired ax */
	delta_ds =s->relvTarget.NearPnt.ds_p- ACC_ECU->DesrDist; /* d_ist-d_soll */
	ax = (delta_ds)/(ACC_ECU->dc_kd)+s->relvTarget.NearPnt.dv_p/ACC_ECU->dc_kv;
	/* ax_sc = desired ax from Speed Control */
	ax_sc =(ACC_ECU->DesrSpd - Vehicle.v)/ ACC_ECU->sc_kv;

	/* Limitation */
	if (ax > ax_sc ) ax = ax_sc;
	if (ax > ACC_ECU->axmax) ax = ACC_ECU->axmax;
	if (ax < ACC_ECU->axmin) ax = ACC_ECU->axmin;
    } else {
	/* Speed Control Algorithm: result = desired ax */
	/*      s->relvTarget.ds = -1; */
	ax = (ACC_ECU->DesrSpd - Vehicle.v) / ACC_ECU->sc_kv;
	/* Limitation */
	if (ax > ACC_ECU->axmax) ax =  ACC_ECU->axmax;
	if (ax < -0.35) ax = -0.35;
    }

    ACC_ECU->DesrAx = ax;

    SetAccelCtrl:
	AccelCtrl.DesrAx = ACC_ECU->DesrAx;
	return;
}


static void
AccelCtrl_Delete (void *MP)
{
    tAccelCtrl *mp = MP;
    memset(mp, 0, sizeof(*mp));
    mp->DesrAx = NOTSET;
}


static void *
AccelCtrl_New (struct tInfos *Inf, const char *kindkey)
{
    char *key, *s, buf[64];
    tACC_ECU *ACC_ECU = &AccelCtrl.ACC_ECU;
    double vInit;

    c_i = 0.0;
    AccelCtrl_Delete (&AccelCtrl);

    /** AccelCtrl */
    key = "AccelCtrl";

    AccelCtrl.p_gain = iGetDblOpt(Inf, strcat(strcpy(buf, key), ".p"), 0.001);
    AccelCtrl.i_gain = iGetDblOpt(Inf, strcat(strcpy(buf, key), ".i"), 1.0);

    /* function pointer for ax-Calculation */
    s = iGetStrOpt(Inf, strcat(strcpy(buf, key), ".DesrAccelFunc"), "ACC");

    if (strcmp(s, "ACC")==0) {
	AccelCtrl.DesrAccelFunc = DesrAccelFunc_ACC;
    } else if (strcmp(s, "DVA")==0) {
	AccelCtrl.DesrAccelFunc = NULL;
    } else if (strcmp(s, "User")==0) {
	if (User_DesrAccelFunc != NULL)
	    AccelCtrl.DesrAccelFunc = User_DesrAccelFunc;
	else
	    AccelCtrl.DesrAccelFunc = NULL;
    } else {
	LogErrF(EC_Init, "AccelCtrl: no supported function '%s' for ax calculation", s);
	return NULL;
    }

    /** ACC */
    key = "AccelCtrl.ACC";

    /* Active / switched on ? */
    ACC_ECU->IsActive = iGetLongOpt (Inf, strcat(strcpy(buf, key), ".IsActive"), 1);

    /* Limit of driver brake to deactivate ACC */
    ACC_ECU->BrakeThreshold = iGetDblOpt(Inf, strcat(strcpy(buf, key), ".BrakeThreshold"), 0.2);

    /* initial time gap / speed */
    ACC_ECU->DesrTGap = iGetDblOpt(Inf, strcat(strcpy(buf, key), ".DesrTGap"), 1.8);
    vInit             = DrivMan.Cfg.Velocity > 10.0*kmh2ms ? DrivMan.Cfg.Velocity : 100*kmh2ms;
    ACC_ECU->DesrSpd  = iGetDblOpt(Inf, strcat(strcpy(buf, key), ".DesrSpd"),  vInit);

    /* controller parameters */
    ACC_ECU->dc_kd = iGetDblOpt(Inf, strcat(strcpy(buf, key), ".DistCtrl.kd"), 36.0);
    ACC_ECU->dc_kv = iGetDblOpt(Inf, strcat(strcpy(buf, key), ".DistCtrl.kv"),  2.0);
    ACC_ECU->sc_kv = iGetDblOpt(Inf, strcat(strcpy(buf, key), ".SpdCtrl.kv"),  13.0);

    /* min/max values */
    ACC_ECU->axmin = iGetDblOpt(Inf, strcat(strcpy(buf, key), ".AxMin"),   -2.5);
    ACC_ECU->axmax = iGetDblOpt(Inf, strcat(strcpy(buf, key), ".AxMax"),    1.0);
    ACC_ECU->dsmin = iGetDblOpt(Inf, strcat(strcpy(buf, key), ".DistMin"), 20.0);
    /* Name of the reference ObjectSensor */
    s = iGetStrOpt(Inf, strcat(strcpy(buf, key), ".RefObjectSensorName"), "RadarL");
    if ((ACC_ECU->RefObjectSensorId = ObjectSensor_FindIndexForName(s))<0) {
	LogErrF(EC_Init, "%s: no ObjectSensor found with the name '%s'", key, s);
	return NULL;
    }

    return &AccelCtrl;
}


static int
AccelCtrl_Calc (void *MP, double dt)
{
    double c, delta_ax, c_p;

    if (SimCore.State != SCState_Simulate || AppStartInfo.ModelCheck ||
	AppStartInfo.DriverAdaption)
	return 0;

    /* Calculate target longitudinal acceleration ax */
    if (AccelCtrl.DesrAccelFunc != NULL)
	AccelCtrl.DesrAccelFunc(dt);

    /* Controller for converting desired ax to gas or brake */
    if (AccelCtrl.DesrAx == NOTSET) {
	/* no control required */
	c_i = VehicleControl.Gas;
	return 0;
    }

    delta_ax = AccelCtrl.DesrAx - Vehicle.PoI_Acc_1[0];
    c_p  = AccelCtrl.p_gain * delta_ax;
    c_i += AccelCtrl.i_gain * delta_ax * dt;
    c = c_p + c_i;	/* PI-Controller */

    /* Limitation */
    if (c >  1) c =  1;
    if (c < -1) c = -1;
    c_i = c - c_p;

    /* Gas or Brake */
    if (c >= 0){
	VehicleControl.Gas   =  c;
	VehicleControl.Brake =  0;
    } else {
	VehicleControl.Gas   =  0;
	VehicleControl.Brake = -c;
    }

    return 0;
}


void
Set_UserDesrAccelFunc(tDesrAccelFunc DesrAccelFunc)
{
    User_DesrAccelFunc = DesrAccelFunc;
}


int
VC_Register_AccelCtrl (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.VehicleControl.New    = 		AccelCtrl_New;
    m.VehicleControl.Calc   = 		AccelCtrl_Calc;
    m.VehicleControl.DeclQuants   = 	AccelCtrl_DeclQuants;
    m.VehicleControl.Delete = 		AccelCtrl_Delete;
    /* Should only be used if the model doesn't read params from extra files */
    m.VehicleControl.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_RegisterIPG(ModelClass_VehicleControl, "AccelCtrl", &m);
}
