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
** Simple 'user vehicle' model of a car
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    Vehicle_Register_MyModel ();
**
******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <math.h>

#include "road.h"
#include "CarMaker.h"
#include "MathUtils.h"
#include "ModelManager.h"
#include "Car/Brake.h"
#include "Vehicle/Tire.h"
#include "Car/PowerTrain.h"

#include "MyModels.h"

#define NPTTRQ 10

static struct {
    double	dt;

    struct {
	double	Dist;
	double	WhlBaseF;
	double	WhlBaseR;
    } Info;

    struct {
	double	Throttle;
	double	Clutch;
	double	Brake;
	int	GearNo;
    } Driver;

    struct {
	tRoadEval *RE;
	int 	LinkObjId, onJunction;
	int 	JuncObjId, nextJuncObjId;
	double	s2nextJunc, s2lastJunc;
	double	Dist;
	double	X_H[3];
	double	Y_H[3];
	double	Z_H[3];
    } Road;

    struct {
	double mass;
	double G;
    } Bdy1;

    struct {
	double x, y, z;

	double vx, vy;
	double vx_1;

	double ax, ay;
	double ax_1, ay_1;

	double Yaw, Yawp, Yawpp;

	double rx, ry, rz;	/* animation */

    } Fr1;

    struct {
	double x, y, z;		/* "driver position" */
    } PoI;

    struct tAero {
	double	rho;
	double	cw;
	double	ca_f;
	double	ca_r;
	double	Ayz;

	double	Fx;
	double	Fz;
    } Aero;

    struct tSteering {
	double i;
	double Ang;
    } Steering;

    struct tWheel {
	double	Radius;
	double	rot;
	double	front_rz;
    } Whl;


    struct tTyre {
	double ax0, ay0;
	double ax_max, ay_max;

	double SideSlipAngF, SideSlipAngR;
	double LongSlipF, LongSlipR;
    } Tyres;

    struct tMyPowerTrain {
	double	rotp[NPTTRQ];
	double	Trq[NPTTRQ];

	double	GearNo;
	double	iAct;
	double	Engine_rotv;
	double	Engine_Trq;
    } PT;
} Car;


/* Skip calculation of the vehicle during preprocessing, just like in CM4SL? */
static int SkipPreprocessing = 0;

/* Vehicle tire data.
   Tires will always be initialized according to the testrun/vehicle dataset,
   but whether they are actually used (i.e. calculated) depends entirely on the
   presence code invoking VhclModel_Tire_Calc() in this module. */
static tTire Tires[VEHICLE_NWHLS];


static void
VehicleModel_DeclQuants (void *MP)
{
}


static void *
VehicleModel_New (struct tInfos *Inf)
{
    struct tVehicleCfg *Cfg = &Vehicle.Cfg;
    int i=0, n, outMask;

    memset (&Car, 0, sizeof(Car));

    /* Get vehicle configuration parameters from vehicle infofile */
    if (Vehicle_GetVhclCfg (Inf, VhclClass_Car_Id) < 0)
	return NULL;

    /* Initialize tires */
    int TireNo, TireFailure = 0;
    for (TireNo=0; TireNo<Vehicle.Cfg.nWheels; TireNo++) {
	const int side = TireNo%2==0 ? 1 : -1;	/* 1 left, -1 right */
	tTire *tire = &Tires[TireNo];

	if (VhclModel_Tire_New(tire, Inf, TireNo, side, 0) != 0)
	    TireFailure = 1;
    }
    if (TireFailure)
	return NULL;

    /* Register Vehicle.Fr1A for sensors and animation. */
    BdyFrame_Register (BdyFrame_MountFrameStr(Vehicle_Fr1A), &Vehicle.Fr1A);

    Car.Info.WhlBaseF = Cfg->CoG2AxleFront;
    Car.Info.WhlBaseR = Cfg->CoG2AxleRear;

    Car.Bdy1.mass =	Cfg->MassTotal;
    Car.Bdy1.G 	=	Car.Bdy1.mass * 9.81;

    Car.Whl.Radius = 	Cfg->WhlRadius;

    Car.Steering.i =	iGetDblOpt(Inf, "Steering.i",	10.0);

    Car.Aero.rho =	iGetDblOpt(Inf, "Aero.rho",	1.205);
    Car.Aero.cw =	iGetDblOpt(Inf, "Aero.cw",	0.91);
    Car.Aero.ca_f =	iGetDblOpt(Inf, "Aero.ca_f",	1.04);
    Car.Aero.ca_r =	iGetDblOpt(Inf, "Aero.ca_r",	1.44);
    Car.Aero.Ayz =	iGetDblOpt(Inf, "Aero.Ayz",	1.725);

    Car.Tyres.ax0 =	iGetDblOpt(Inf, "Tire.ax0",	13.0);
    Car.Tyres.ay0 =	iGetDblOpt(Inf, "Tire.ay0",	13.0);

    n = 0;
    double *dv = iGetTableOpt2 (Inf, "PowerTrain.Engine.Trq", NULL, 2, &n);
    if (dv != NULL && n > 0) {
	n = M_MIN(NPTTRQ, n);
	for (i=0; i<n; i++) {
	    Car.PT.rotp[i] = dv[i] * rpm2radsec;
	    Car.PT.Trq[i]  = dv[i+n];
	}
	free (dv);
    }

    for (; i < NPTTRQ; i++) {
	Car.PT.rotp[i] =   Car.PT.rotp[i-1] + 1*rpm2radsec;
	Car.PT.Trq[i] =	   Car.PT.Trq[i-1];
    }

    /* Initial conditions. */
    Car.Fr1.vx_1  = DrivMan.Cfg.Velocity;
    Car.Road.Dist = Vehicle.sRoad = DrivMan.Cfg.sRoad;
    Car.Fr1.x	  = Cfg->Fr1_Pos[0];
    Car.Fr1.y	  = Cfg->Fr1_Pos[1];
    Car.Fr1.z	  = Cfg->Fr1_Pos[2];
    Car.Fr1.Yaw	  = Cfg->Yaw;
    Car.Fr1.vx    = Car.Fr1.vx_1 * cos(Car.Fr1.Yaw);
    Car.Fr1.vy    = Car.Fr1.vx_1 * sin(Car.Fr1.Yaw);


    /* Provide vehicle geometry information for animation tool IPGMovie. */
    if (Vehicle_CreateVhclMsgFromCfg (VhclClass_Car_Id) < 0)
	return NULL;

    /* RoadEval Handle (only on valid route) */
    if (Env.Route.ObjId < 0) {
	LogErrF (EC_Init, "MyCar: disabled (works only on route)");
	return NULL;
    }

    outMask = ROAD_OT_MIN | ROAD_OT_SUV | ROAD_OT_TUV | ROAD_OT_NUV | ROAD_OT_EXT;
    if ((Car.Road.RE = RoadNewRoadEval (Env.Road, ROAD_BUMP_CRG, outMask, "MyCar")) == NULL) {
	LogErrF (EC_Init, "MyCar: failed to get road handle");
	return NULL;
    }
    RoadEvalSetRouteByObjId (Car.Road.RE, Env.Route.ObjId, 0);


#if 0
    SkipPreprocessing = iGetBoolOpt(Inf, "Vehicle.SkipPreprocessing", 0);
#else
    // The MyCar.c vehicle model must never be invoked during pre- or postprocessing.
    SkipPreprocessing = 1;
#endif

    /* Disable postprocessing after end of simulation, i.e. DON'T pilot the
       vehicle to a stand-still state, freeze actual vehicle state instead. */
    SimCore.GetIdle.Skip = 1;

    if (SkipPreprocessing) {
	/* Provide roughly meaningful values for the first invocation of IPGDriver
	   in state "Simulate", as the first call to our vehicle model's calc()-
	   function only comes afterwards when not taking part in preprocessing. */
	Vehicle.PoI_Pos[0] = DrivMan.Cfg.RoadPos[0];
	Vehicle.PoI_Pos[1] = DrivMan.Cfg.RoadPos[1];
	Vehicle.PoI_Pos[2] = DrivMan.Cfg.RoadPos[2];
	Vehicle.PoI_Vel[0] = DrivMan.Cfg.Velocity * DrivMan.Cfg.RoadDir[0];
	Vehicle.PoI_Vel[1] = DrivMan.Cfg.Velocity * DrivMan.Cfg.RoadDir[1];
    }

    /* Our vehicle data (global variable Car) is statically allocated, so we don't
       need a handle to it, but returning NULL would mean failure, so we return a
       dummy address instead. */
    return (void *)1;
}


static void
Car2Vehicle (void)
{
    double val;

    Vehicle.sRoad =		Car.Road.Dist;
    Vehicle.Distance =		Car.Info.Dist;
    Vehicle.v =			Car.Fr1.vx_1;

    /* Information of location on road network (to be considered by autonomous traffic) */
    Vehicle.Road.LinkObjId =	Car.Road.LinkObjId;
    Vehicle.Road.onJunction =	Car.Road.onJunction;
    Vehicle.Road.JuncObjId =	Car.Road.JuncObjId;
    Vehicle.Road.nextJuncObjId = Car.Road.nextJuncObjId;
    Vehicle.Road.s2nextJunc = 	Car.Road.s2nextJunc;
    Vehicle.Road.s2lastJunc =	Car.Road.s2lastJunc;

    Vehicle.PoI_Pos[0] =	Car.PoI.x;
    Vehicle.PoI_Pos[1] =	Car.PoI.y;
    Vehicle.PoI_Pos[2] =	Car.PoI.z;

    Vehicle.PoI_Vel[0] = 	Car.Fr1.vx;
    Vehicle.PoI_Vel[1] = 	Car.Fr1.vy;
    Vehicle.PoI_Vel[2] = 	0.0;

    Vehicle.PoI_Acc[0] =	Car.Fr1.ax;
    Vehicle.PoI_Acc[1] =	Car.Fr1.ay;
    Vehicle.PoI_Acc[2] =	0.0;

    Vehicle.PoI_Vel_1[0] = 	Car.Fr1.vx_1;
    Vehicle.PoI_Vel_1[1] = 	0.0;
    Vehicle.PoI_Vel_1[2] = 	0.0;

    Vehicle.PoI_Acc_1[0] =	Car.Fr1.ax_1;
    Vehicle.PoI_Acc_1[1] =	Car.Fr1.ay_1;
    Vehicle.PoI_Acc_1[2] =	0.0;

    Vehicle.Yaw =		Car.Fr1.Yaw;
    Vehicle.YawRate =		Car.Fr1.Yawp;
    Vehicle.YawAcc =		Car.Fr1.Yawpp;

    Vehicle.Roll =		Car.Fr1.rx;
    Vehicle.Pitch =		Car.Fr1.ry;

    Vehicle.FL.r_zxy[2] =	Car.Whl.front_rz;
    Vehicle.FR.r_zxy[2] =	Car.Whl.front_rz;

    Vehicle.FL.rot =		Car.Whl.rot;
    Vehicle.FR.rot =		Car.Whl.rot;
    Vehicle.RL.rot =		Car.Whl.rot;
    Vehicle.RR.rot =		Car.Whl.rot;

    Vehicle.FL.LongSlip =	Car.Tyres.LongSlipF;
    Vehicle.FR.LongSlip =	Car.Tyres.LongSlipF;
    Vehicle.RL.LongSlip =	Car.Tyres.LongSlipR;
    Vehicle.RR.LongSlip =	Car.Tyres.LongSlipR;

    Vehicle.FL.SideSlip =	Car.Tyres.SideSlipAngF;
    Vehicle.FR.SideSlip =	Car.Tyres.SideSlipAngF;
    Vehicle.RL.SideSlip =	Car.Tyres.SideSlipAngR;
    Vehicle.RR.SideSlip =	Car.Tyres.SideSlipAngR;

    val =			Car.Bdy1.G*0.25;
    Vehicle.FL.Fz =		val;
    Vehicle.FR.Fz =		val;
    Vehicle.RL.Fz =		val;
    Vehicle.RR.Fz =		val;

    Vehicle.Steering.Ang =	Car.Steering.Ang;
    Vehicle.Steering.AngVel =	0.0;
    Vehicle.Steering.AngAcc =	0.0;

    /* Vehicle Fr1A for sensors and animation */
    Vehicle.Fr1A.t_0[0] = Car.Fr1.x;
    Vehicle.Fr1A.t_0[1] = Car.Fr1.y;
    Vehicle.Fr1A.t_0[2] = Car.Fr1.z;

    Vehicle.Fr1A.v_0[0] = Car.Fr1.vx;
    Vehicle.Fr1A.v_0[1] = Car.Fr1.vy;
    Vehicle.Fr1A.v_0[2] = 0.0;

    Vehicle.Fr1A.a_0[0] = Car.Fr1.ax;
    Vehicle.Fr1A.a_0[1] = Car.Fr1.ay;
    Vehicle.Fr1A.a_0[2] = 0.0;

    if (GCSOn) {
	GCS_ConvFr0toGCS(Vehicle.PoI_Pos, &Vehicle.PoI_GCS);
    } else {
	Vehicle.PoI_GCS.Elev = Vehicle.PoI_Pos[2];
    }

#define Calc_Tr2Fr0_ZYX(Tr2FrA, rx, ry, rz) \
    do { \
	    double srx, crx, sry, cry, srz, crz; \
	    \
	    M_SINCOS((rx), &srx, &crx); \
	    M_SINCOS((ry), &sry, &cry); \
	    M_SINCOS((rz), &srz, &crz); \
	    \
	    (Tr2FrA)[0][0] =  cry*crz; \
	    (Tr2FrA)[0][1] =  srx*sry*crz - crx*srz; \
	    (Tr2FrA)[0][2] =  crx*sry*crz + srx*srz; \
	    (Tr2FrA)[1][0] =  cry*srz; \
	    (Tr2FrA)[1][1] =  srx*sry*srz + crx*crz; \
	    (Tr2FrA)[1][2] =  crx*sry*srz - srx*crz; \
	    (Tr2FrA)[2][0] = -sry; \
	    (Tr2FrA)[2][1] =  srx*cry; \
	    (Tr2FrA)[2][2] =  crx*cry; \
    } while (0)

    Calc_Tr2Fr0_ZYX (Vehicle.Fr1A.Tr2Fr0, Car.Fr1.rx, Car.Fr1.ry, Car.Fr1.Yaw);

    VEC_AssignFromMatCol (Vehicle.Fr1A.X_0, Vehicle.Fr1A.Tr2Fr0, 0);
    VEC_AssignFromMatCol (Vehicle.Fr1A.Y_0, Vehicle.Fr1A.Tr2Fr0, 1);
    VEC_AssignFromMatCol (Vehicle.Fr1A.Z_0, Vehicle.Fr1A.Tr2Fr0, 2);

    PowerTrain.IF.OperationState =	OperState_Driving;
    PowerTrain.IF.Engine_rotv =		Car.PT.Engine_rotv;
    PowerTrain.IF.GearNo =		Car.Driver.GearNo;

    PowerTrain.IF.WheelOut[0].rot =	Car.Whl.rot;
    PowerTrain.IF.WheelOut[1].rot =	Car.Whl.rot;
    PowerTrain.IF.WheelOut[2].rot =	Car.Whl.rot;
    PowerTrain.IF.WheelOut[3].rot =	Car.Whl.rot;
}


static int
VehicleModel_Calc (void *MP, double dt)
{;
    double val;
    const char *msgpre = "MyCar: ";
    tRoadRouteIn     rIn;
    tRoadRouteOutExt rOut;
    struct tVehicleCfg *Cfg = &Vehicle.Cfg;

    Car.dt = dt;

#if 0
    if (SkipPreprocessing) {
	if (SCState_Start<=SimCore.State && SimCore.State<=SCState_StartLastCycle)
	    goto OkReturn;
    }
#else
    // The MyCar.c vehicle model must never be invoked during pre- or postprocessing.
    if (SimCore.State != SCState_Simulate)
	goto OkReturn;
#endif

    /*** driver activity to vehicle */
    Car.Driver.GearNo =		VehicleControl.GearNo;
    Car.Driver.Throttle =	VehicleControl.Gas;
    Car.Driver.Clutch =		VehicleControl.Clutch;
    Car.Driver.Brake =		VehicleControl.Brake;
    Car.Steering.Ang =		VehicleControl.Steering.Ang;

    /*** vehicle on the road */
    memset(&rIn, 0, sizeof(struct tRoadRouteIn));
    VEC_AssignByComp (rIn.xyz, Car.Fr1.x, Car.Fr1.y, Car.Fr1.z);
    rIn.st[0] = Car.Road.Dist;

    if (RoadRouteEvalExt (Car.Road.RE, msgpre, RIT_XY_S, &rIn, &rOut) != ROAD_Ok)
	goto ErrorReturn;
    if (!rOut.onRoad) {
	LogErrF (EC_Sim, "MyCar leaves road at about sRoad=%g m, x=%g, y=%g",
		 Car.Road.Dist, rOut.xyz[0], rOut.xyz[1]);
	goto ErrorReturn;
    }

    Car.Fr1.z     		= rOut.xyz[2];
    Car.Road.Dist 		= rOut.st[0];
    Car.Road.LinkObjId 		= rOut.linkObjId;
    Car.Road.onJunction 	= rOut.onJunction;
    Car.Road.JuncObjId 		= rOut.juncObjId;
    Car.Road.nextJuncObjId 	= rOut.nextJuncObjId;
    Car.Road.s2nextJunc 	= rOut.nextJuncDist;
    Car.Road.s2lastJunc 	= rOut.lastJuncDist;

    VEC_Assign (Car.Road.X_H, rOut.suv);
    VEC_Assign (Car.Road.Y_H, rOut.tuv);
    VEC_Assign (Car.Road.Z_H, rOut.nuv);


    /*** engine torque */
    {
	int i, GearNo = M_MAX(Car.Driver.GearNo, 1);
	double dTrq_dn, dn;

	Car.PT.iAct = Cfg->iDiff * Cfg->iFGear[GearNo];
	if (Car.PT.iAct != 0.0) {
	    Car.PT.Engine_rotv = Car.Fr1.vx_1 * Car.PT.iAct / Car.Whl.Radius;
	}
	if (Car.PT.Engine_rotv < Car.PT.rotp[0])
	    Car.PT.Engine_rotv = Car.PT.rotp[0];

	for (i=1; i < NPTTRQ; i++) {
	    if (Car.PT.Engine_rotv < Car.PT.rotp[i])
		break;
	}

	i = M_MIN(NPTTRQ-1, i);
	dn      = Car.PT.Engine_rotv - Car.PT.rotp[i-1];
	dTrq_dn = (Car.PT.Trq[i]  - Car.PT.Trq[i-1])
		/ (Car.PT.rotp[i] - Car.PT.rotp[i-1]);

	Car.PT.Engine_Trq = Car.Driver.Throttle * (Car.PT.Trq[i-1] + dTrq_dn*dn);
    }

    /*** aerodynamics */
    val = 0.5 * Car.Aero.rho * Car.Aero.Ayz * (Car.Fr1.vx_1*Car.Fr1.vx_1);
    Car.Aero.Fx = val * Car.Aero.cw;
    Car.Aero.Fz = val * (Car.Aero.ca_f+Car.Aero.ca_r);


#if 0
    int TireNo;

    /* This is just a rough sketch of how to use the built-in tire models,
       as this particular vehicle model provides its own tire calculation. */

    for (TireNo=0; TireNo<Vehicle.Cfg.nWheels; TireNo++) {
	tTire *tire = &Tires[TireNo];
	tTireIF *IF = &tire->IF;

	/* MISSING HERE: Set tire inputs */
	if (tire->md.Tire.is3DTire) {
	    // IF->WC_t_0[0] = ...;
	    // IF->WC_t_0[1] = ...;
	    // IF->WC_t_0[2] = ...;
	    // ...
	} else {
	    // IF->Frc_W[0] = ...;
	    // IF->Frc_W[1] = ...;
	    // IF->Frc_W[2] = ...;
	    // ...
	}

	VhclModel_Tire_Calc(SimCore.DeltaT, tire);

	/* MISSING HERE: Fetch tire outputs */
	if (tire->md.Tire.is3DTire) {
	    // ... = IF->Frc_C[0];
	    // ... = IF->Frc_C[1];
	    // ... = IF->Frc_C[2];
	    // ...
	} else {
	    // ... = IF->Slp;
	    // ... = IF->Alpha;
	    // ... = IF->TurnSlp;
	    // ...
	}
    }
#else
    /*** tyre forces */
    Car.Tyres.ax_max = Car.Tyres.ax0 * (1. + Car.Aero.Fz/Car.Bdy1.G);
    Car.Tyres.ay_max = Car.Tyres.ay0 * (1. + Car.Aero.Fz/Car.Bdy1.G);

    {
	const double vCritical = 2000.0 / 3.6;
	double facBrake = Car.Fr1.vx_1<=0 ? 0.0 : 1.0;
	Car.Fr1.ax_1
	    = Car.PT.Engine_Trq * Car.PT.iAct / Car.Whl.Radius / Car.Bdy1.mass
	    - 40. * Car.Driver.Brake * facBrake;

	Car.Fr1.ay_1
	    = Car.Fr1.vx_1*Car.Fr1.vx_1
	    * Car.Steering.Ang/Car.Steering.i
	    * (1.0/(Car.Info.WhlBaseF-Car.Info.WhlBaseR)
	       - Car.Fr1.vx_1/vCritical);
    }

    /*** longitudinal and side slip */
    val = 1.2*Car.Fr1.ax_1/Car.Tyres.ax_max;
    if      (val >  1.) Car.Tyres.LongSlipF =  0.05 + 10. * (val-1.0);
    else if (val < -1.) Car.Tyres.LongSlipF = -0.05 + 10. * (val+1.0);
    else                  Car.Tyres.LongSlipF =  0.05 * asin(val);

    if      (Car.Tyres.LongSlipF >   1.) Car.Tyres.LongSlipF =  1.0;
    else if (Car.Tyres.LongSlipF <  -1.) Car.Tyres.LongSlipF = -1.0;
    Car.Tyres.LongSlipR = Car.Tyres.LongSlipF;


    val = 1.2*Car.Fr1.ay_1/Car.Tyres.ay_max;
    if      (val >  1.) Car.Tyres.SideSlipAngF =  0.05 + 10. * (val-1.0);
    else if (val < -1.) Car.Tyres.SideSlipAngF = -0.05 + 10. * (val+1.0);
    else                Car.Tyres.SideSlipAngF =  0.05 * asin(val);

    if      (Car.Tyres.SideSlipAngF >  M_PI/2.) Car.Tyres.SideSlipAngF =  M_PI/2.;
    else if (Car.Tyres.SideSlipAngF < -M_PI/2.) Car.Tyres.SideSlipAngF = -M_PI/2.;
    Car.Tyres.SideSlipAngR = Car.Tyres.SideSlipAngF;

    /* tyre force limitation */
    if (Car.Fr1.ay_1 >  Car.Tyres.ay_max) Car.Fr1.ay_1 =  Car.Tyres.ay_max;
    if (Car.Fr1.ay_1 < -Car.Tyres.ay_max) Car.Fr1.ay_1 = -Car.Tyres.ay_max;

    if (Car.Fr1.ax_1 >  Car.Tyres.ax_max) Car.Fr1.ax_1 =  Car.Tyres.ax_max;
    if (Car.Fr1.ax_1 < -Car.Tyres.ax_max) Car.Fr1.ax_1 = -Car.Tyres.ax_max;
#endif


    /* Aerodynamik */
    /* =========== */
    /* Hier erst, da nicht ueber Reifen uebertragen */
    Car.Fr1.ax_1 -= Car.Aero.Fx / Car.Bdy1.mass;


    /* vehicle motion */
    Car.Fr1.vx_1 += Car.Fr1.ax_1*dt;


    Car.Fr1.ax
	= Car.Fr1.ax_1 * cos(Car.Fr1.Yaw) - Car.Fr1.ay_1 * sin(Car.Fr1.Yaw);
    Car.Fr1.ay
	= Car.Fr1.ax_1 * sin(Car.Fr1.Yaw) + Car.Fr1.ay_1 * cos(Car.Fr1.Yaw);

    Car.Fr1.vx += Car.Fr1.ax * dt;
    Car.Fr1.vy += Car.Fr1.ay * dt;

    Car.Fr1.x += Car.Fr1.vx * dt + 0.5 * Car.Fr1.ax * dt * dt;
    Car.Fr1.y += Car.Fr1.vy * dt + 0.5 * Car.Fr1.ay * dt * dt;

    Car.Fr1.Yawpp =	0.0;
    if (Car.Fr1.vx_1 > 0.01 || Car.Fr1.vx_1 < -0.01)
	Car.Fr1.Yawp =	Car.Fr1.ay_1/Car.Fr1.vx_1;
    else
	Car.Fr1.Yawp =	0.0;
    Car.Fr1.Yaw +=	Car.Fr1.Yawp * dt;

    {
	/* Calculating the PoI position (driver point in CoM) */
	double Pos01PoI_0[2];

	/* vector Fr1->PoI, expressed in Fr0 */
	Pos01PoI_0[0] = Vehicle.Cfg.Bdy1_CoM[0] * cos(Car.Fr1.Yaw) -
			Vehicle.Cfg.Bdy1_CoM[1] * sin(Car.Fr1.Yaw);
	Pos01PoI_0[1] = Vehicle.Cfg.Bdy1_CoM[0] * sin(Car.Fr1.Yaw) +
			Vehicle.Cfg.Bdy1_CoM[1] * cos(Car.Fr1.Yaw);

	/* global PoI position, expressed in Fr0 */
	Car.PoI.x = Car.Fr1.x + Pos01PoI_0[0];
	Car.PoI.y = Car.Fr1.y + Pos01PoI_0[1];
	Car.PoI.z = Car.Fr1.z;
    }

    /*** calculate additional output signals	*/

    /*** milometer */
    Car.Info.Dist +=	Car.Fr1.vx_1 * Car.dt;

    /*** animation quantities: joint angles, used by IPGMovie */
    NEV2FreiZYX (Car.Road.X_H, Car.Road.Y_H, Car.Road.Z_H,
		 &Car.Fr1.rx, &Car.Fr1.ry, &val);
    Car.Fr1.rz =	Car.Fr1.Yaw;

    Car.Whl.rot +=	Car.Fr1.vx_1 / (2.0 * M_PI * Car.Whl.Radius) * Car.dt;
    Car.Whl.front_rz =	Car.Steering.Ang / Car.Steering.i;

    /*** assign values from Car to the Vehicle struct */
    Car2Vehicle ();

  OkReturn:
    return 0;

  ErrorReturn:
    return -1;
}


static void
VehicleModel_Delete (void *MP)
{
    int TireNo;
    for (TireNo=0; TireNo<Vehicle.Cfg.nWheels; TireNo++)
	VhclModel_Tire_Delete(TireNo);

    if (Car.Road.RE != NULL) {
	RoadDeleteRoadEval (Car.Road.RE);
	Car.Road.RE = NULL;
    }
}


int
Vehicle_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.Vehicle.VersionId =		1;
    m.Vehicle.CompatVersionId =		1;
    m.Vehicle.New =			VehicleModel_New;
    m.Vehicle.Calc =			VehicleModel_Calc;
    m.Vehicle.Delete =			VehicleModel_Delete;
    m.Vehicle.DeclQuants =		VehicleModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.Vehicle.ParamsChanged = 		ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_Vehicle, "MyCar", &m);
}
