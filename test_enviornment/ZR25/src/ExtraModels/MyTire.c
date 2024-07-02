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
** (turning slip and inclination angle are not considered)
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    Tire_Register_MyModelCPI ();
**    Tire_Register_MyModelSTI ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "Tire";

/* User defined model struct */
struct tMyModel {
    /* Parameters */
    tTParamHeader Head;		/* Struct MUST begin with the tire header struct */

    struct {			/* Fx, Fy forces curve coefficients */
	double	c1;
	double	c2;
	double	c3;
    } Fx, Fy;

    double	cMy;		/* Rolling resistance torque coefficient */

    double	defl;		/* Tire deflection */
    double	rlen;		/* Tire relaxation length [m] */
    double	Slip;		/* Long. Slip [] */
    double	Alpha;		/* Side Slip [rad] */
};


static void
MyModel_DeclQuants_dyn (struct tMyModel *mp, int park)
{
    static struct tMyModel MyModel_Dummy = {{0}};
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


static void
MyModel_Delete (void *MP)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Park the dict entries for dynamically allocated variables before deleting */
    MyModel_DeclQuants_dyn (mp, 1);

    if (mp->Head.ttparam != NULL)
	TireTool_Delete (mp->Head.ttparam);

    free(mp);
}



/* Tire for Contact Point Interface (CPI) ************************************/

static const char ThisModelKind_CPI[] = "MyModelCPI";
static const int  ThisVersionId_CPI   = 1;

static void *
MyModel_CPI_New (
    struct tInfos	*inf,
    const char		*ParamDir,
    int			TireNo,
    int			Side,
    int			Twin)
{
    struct tMyModel *mp = NULL;
    char MsgPre[64];

    mp = (struct tMyModel*)calloc(1,sizeof(*mp));

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind_CPI);

    mp->Fx.c1 = iGetDbl (inf, "MyTire.Fx.c1");
    mp->Fx.c2 = iGetDbl (inf, "MyTire.Fx.c2");
    mp->Fx.c3 = iGetDbl (inf, "MyTire.Fx.c3");

    mp->Fy.c1 = iGetDbl (inf, "MyTire.Fy.c1");
    mp->Fy.c2 = iGetDbl (inf, "MyTire.Fy.c2");
    mp->Fy.c3 = iGetDbl (inf, "MyTire.Fy.c3");

    mp->cMy =   iGetDbl (inf, "MyTire.cMy");
    mp->rlen =  iGetDbl (inf, "MyTire.rlen");

    /* Tire header initialization Tire_Head_Init()
       and Tire_New_Head MUST be called */
    Tire_Head_Init (&mp->Head);

    if (Tire_New_Head(&mp->Head, inf, TireNo, Side, Twin) < 0) {
	LogErrF (EC_Init, "%s: Failed to initialize the tire header struct", MsgPre);
	goto ErrorReturn;
    }

    /* Get parameters for the standstill model */
    if ((mp->Head.ttparam = TireTool_New(inf)) == NULL) {
	LogErrF (EC_Init, "%s: Failed to initialize the standstill model", MsgPre);
	goto ErrorReturn;
    }

    return mp;

     ErrorReturn:
	free (mp);
	return NULL;
}


static void
Calc_TireFrcTrq_W (struct tMyModel *mp, tTireIF *IF, double dt)
{
    double T, abss, absa;

    /* Relaxation of the slips */
    if (dt > 0) {
	T = mp->rlen / M_MAX(fabs(IF->P_v0_W[0]), 2.0);
	FILTER_LP2 (IF->Slp,   mp->Slip,  T, dt);
	FILTER_LP2 (IF->Alpha, mp->Alpha, T, dt);
    }

    /* Fx */
    abss = fabs(IF->Slp);
    IF->Frc_W[0]  = IF->Frc_W[2] * M_SGN(IF->Slp) * mp->Fx.c1 * 
    		    (1.0 - exp(-mp->Fx.c2*abss)) - mp->Fx.c3*abss;
    IF->Frc_W[0] *= IF->muRoad;

    /* Fy */
    absa = fabs(sin(IF->Alpha));
    IF->Frc_W[1]  = -IF->Frc_W[2] * M_SGN(IF->Alpha) * mp->Fy.c1 * 
    		    (1.0 - exp(-mp->Fy.c2*absa)) - mp->Fy.c3*absa;
    IF->Frc_W[1] *= IF->muRoad;

    /* My */
    IF->Trq_W[1]  = IF->Frc_W[2] * mp->cMy * IF->Radius * (IF->vBelt>0 ? -1.0 : 1.0);
    IF->Trq_W[1] *= IF->muRoad;
}


static int
MyModel_CPI_Calc (void *MP, tTireIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    if (IF->TireCalcMode == TireCalcPrepNoRoad ||
	IF->TireCalcMode == TireCalcPrepWithRoad) {
	//override the IF->Frc_W[2] for preparation phase
	return 0;
    }

    /* Calculate the slip and the side slip with the standstill model */
    TireTool_Kin_Calc (mp->Head.ttparam, IF, dt);

    /* Kinematic rolling radius */
    IF->rBelt_eff = mp->Head.KinRollRadius;

    /* Tire forces and torque in Frame W */
    Calc_TireFrcTrq_W (mp, IF, dt);

    return 0;
}


int 
Tire_Register_MyModelCPI (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.Tire.VersionId =		ThisVersionId_CPI;
    m.Tire.is3DTire = 		0;	// set 0 if CPI interface is used
    m.Tire.New =		MyModel_CPI_New;
    m.Tire.Calc =		MyModel_CPI_Calc;
    m.Tire.Delete =		MyModel_Delete;
    m.Tire.DeclQuants =		MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.Tire.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_Tire, ThisModelKind_CPI, &m);
}


/* Tire for Standard Tire Interface (STI) ************************************/

static const char ThisModelKind_STI[] = "MyModelSTI";
static const int  ThisVersionId_STI   = 1;

static void *
MyModel_STI_New (
    struct tInfos	*inf,
    const char		*ParamDir,
    int			TireNo,
    int			Side,
    int			Twin)
{
    struct tMyModel *mp = NULL;
    char MsgPre[64];

    mp = (struct tMyModel*)calloc(1,sizeof(*mp));

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind_STI);

    mp->Fx.c1 = iGetDbl (inf, "MyTire.Fx.c1");
    mp->Fx.c2 = iGetDbl (inf, "MyTire.Fx.c2");
    mp->Fx.c3 = iGetDbl (inf, "MyTire.Fx.c3");

    mp->Fy.c1 = iGetDbl (inf, "MyTire.Fy.c1");
    mp->Fy.c2 = iGetDbl (inf, "MyTire.Fy.c2");
    mp->Fy.c3 = iGetDbl (inf, "MyTire.Fy.c3");

    mp->cMy =   iGetDbl (inf, "MyTire.cMy");
    mp->rlen =  iGetDbl (inf, "MyTire.rlen");

    /* Tire header initialization Tire_Head_Init()
       and Tire_New_Head MUST be called */
    Tire_Head_Init (&mp->Head);

    if (Tire_New_Head(&mp->Head, inf, TireNo, Side, Twin) < 0) {
	LogErrF (EC_Init, "%s: Failed to initialize the tire header struct", MsgPre);
	goto ErrorReturn;
    }

    /* Get parameters for the standstill model */
    if ((mp->Head.ttparam = TireTool_New(inf)) == NULL) {
	LogErrF (EC_Init, "%s: Failed to initialize the standstill model", MsgPre);
	goto ErrorReturn;
    }

    return mp;

  ErrorReturn:
	free (mp);
	return NULL;
}

static void
Set_RoadOut_ForNoRoad (tRoadGeoIn *rIn, tRoadGeoOut *rOut)
{
    rOut->xyz[0]  = rIn->xyz[0];
    rOut->xyz[1]  = rIn->xyz[1];
    VEC_Assign (rOut->nuv, EZ3x1);
    rOut->fric    = 1.0;
}


static void
MyModel_STI_Calc_CP (tTParamHeader *th, tTireIF *IF)
{
    int    i, rv = 0;
    double Y_C[3], X_W[3], Y_W[3], Z_W[3], CR[3], CP[3], R_0[3];
    double alpha, beta, wxr_0[3], vP_0[3], omega_0[3], vel_0[3];
    tRoadGeoIn  rIn;
    tRoadGeoOut rOut;

    /* Wheel spin axis */
    Y_C[0] = IF->WC_Tr2Fr0[0][1];
    Y_C[1] = IF->WC_Tr2Fr0[1][1];
    Y_C[2] = IF->WC_Tr2Fr0[2][1];

    /* Contact point estimation: new WC_t_0 + last CP_0 */
    memset(&rIn, 0, sizeof(struct tRoadGeoIn));
    VEC_Add (rIn.xyz, IF->WC_t_0, IF->CP_0);

    if ((IF->TireCalcMode == TireCalcPrepWithRoad
      || IF->TireCalcMode == TireCalcPrepWhlWithRoad
      || IF->TireCalcMode == TireCalcSim) && !Vehicle.ModelCheck.SPMM) {

	rv = RoadGeoEval (th->RoadEval, NULL, &rIn, &rOut);
	if (rv != ROAD_Ok) {
	    rOut.xyz[2] = IF->P_0[2]; /* last z-Value */
	    rOut.onRoad = RPK_OffRoad;
	    Set_RoadOut_ForNoRoad (&rIn, &rOut);
	    rv = 1;
	}

    } else {
	/* Flat plane road for preparation phase */
	rOut.xyz[2] = 0.0;
	rOut.onRoad = RPK_OnRoad;
	Set_RoadOut_ForNoRoad (&rIn, &rOut);
    }

    VEC_Assign (R_0, rOut.xyz);
    IF->muRoad = rOut.fric;
    IF->onRoad = rOut.onRoad;

    VEC_Assign (Z_W, rOut.nuv);

    /* Tire frame W */
    VEC_Cross     (X_W, Y_C, Z_W);
    VEC_Normalize (X_W, X_W);
    VEC_Cross     (Y_W, Z_W, X_W);

    IF->FrW_Tr2Fr0[0][0] = X_W[0];
    IF->FrW_Tr2Fr0[1][0] = X_W[1];
    IF->FrW_Tr2Fr0[2][0] = X_W[2];

    IF->FrW_Tr2Fr0[0][1] = Y_W[0];
    IF->FrW_Tr2Fr0[1][1] = Y_W[1];
    IF->FrW_Tr2Fr0[2][1] = Y_W[2];

    IF->FrW_Tr2Fr0[0][2] = Z_W[0];
    IF->FrW_Tr2Fr0[1][2] = Z_W[1];
    IF->FrW_Tr2Fr0[2][2] = Z_W[2];

    /* Contact point correction */
    VEC_Sub (CR, R_0, IF->WC_t_0);
    alpha = - VEC_Scalar(CR, X_W);
    beta  = - VEC_Scalar(CR, Y_C) / VEC_Scalar(Y_W, Y_C);
    for (i=0; i < 3; i++)
	CP[i] = CR[i] + alpha * X_W[i] + beta * Y_W[i];

    if (IF->onRoad == RPK_OffRoad && rv >= 0) {
	LogErrF (EC_Sim, "Vehicle leaves road at x=%g, y=%g TireNo=%d",
	                 rIn.xyz[0], rIn.xyz[1], th->TireNo);
    }

    VEC_MatVec (vel_0,   IF->WC_Tr2Fr0, IF->WC_vel_2);
    VEC_MatVec (omega_0, IF->WC_Tr2Fr0, IF->WC_omega_2);

    VEC_Cross   (wxr_0, omega_0, CP);
    VEC_Add     (vP_0, vel_0, wxr_0);
    VEC_MatTVec (IF->P_v0_W, IF->FrW_Tr2Fr0, vP_0);

    VEC_Add    (IF->P_0, IF->WC_t_0, CP);
    VEC_Assign (IF->CP_0,  CP);
    VEC_Assign (IF->Z_W_0, Z_W);
}


static void
MyModel_STI_Calc_Frc2WC (tTireIF *IF)
{
    double v0[3], v1[3], v2[3], WCP_0[3];

    /* Vector wheel carrier center to contact point considering
       lateral offset for twin tires */
    VEC_Add (WCP_0, IF->WC2Tire_0, IF->CP_0);

    /* Forces from contact point to wheel carrier center */
    VEC_MatVec	(v0, 		IF->FrW_Tr2Fr0,	IF->Frc_W);
    VEC_MatTVec	(IF->Frc_C,	IF->WC_Tr2Fr0,	v0);

    /* Torque from contact point to wheel carrier center */
    VEC_Cross	(v2, WCP_0, v0);

    VEC_MatVec	(v1, 		IF->FrW_Tr2Fr0,	IF->Trq_W);
    VEC_Add	(v1, v1, v2);
    VEC_MatTVec	(IF->Trq_C,	IF->WC_Tr2Fr0,	v1);
}


static int
MyModel_STI_Calc (void *MP, tTireIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    tTParamHeader *th   = &mp->Head;
    int i;
    double defl_last, deflv;

    /* Modelcheck: Fz and muRoad comes as input into the model */
    if (IF->TireCalcMode == TireCalcModelCheck)
	goto FrcTrq;

    /* Set tire forces and torques to zero */
    for (i=0; i<3; i++) {
	IF->Frc_W[i] = 0.0;
	IF->Trq_W[i] = 0.0;
    }

    /** Contact point calculation */
    MyModel_STI_Calc_CP (th, IF);

    /* Loaded tire radius */
    IF->Radius = VEC_Norm (IF->CP_0);

    /* Normal force Fz */
    defl_last = mp->defl;
    mp->defl  = th->NomRadius + M_COPYSIGN(IF->Radius, VEC_Scalar(IF->Z_W_0, IF->CP_0));
    deflv = (dt > 0) ? (mp->defl - defl_last) / dt : 0.0;
    IF->Frc_W[2] = mp->defl * th->Radial.Stiffness - deflv * th->Radial.Damping;

    if (IF->TireCalcMode == TireCalcPrepNoRoad ||
	IF->TireCalcMode == TireCalcPrepWithRoad) {
	//only Fz calculation while preparation phase
	goto Frc2WC;
    }

    if (IF->Frc_W[2] < 0) {
	//no forces and torques with negative normal force
	goto Frc2WC;
    }

    /** Tire forces and torque in contact point */
    FrcTrq:
    /* Kinematic rolling radius */
    IF->rBelt_eff = mp->Head.KinRollRadius;
    IF->vBelt     = IF->rBelt_eff * IF->Rim_rotv;

    /* Calculate the slip and the side slip with the standstill model */
    TireTool_Kin_Calc (mp->Head.ttparam, IF, dt);

    /* Tire forces and torque in Frame W */
    Calc_TireFrcTrq_W (mp, IF, dt);


    /** Tire forces and torques in wheel carrier center */
    Frc2WC:
    MyModel_STI_Calc_Frc2WC (IF);

    return 0;
}


int 
Tire_Register_MyModelSTI (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.Tire.VersionId =		ThisVersionId_STI;
    m.Tire.is3DTire = 		1;	// set 1 if STI interface is used
    m.Tire.New =		MyModel_STI_New;
    m.Tire.Calc =		MyModel_STI_Calc;
    m.Tire.Delete =		MyModel_Delete;
    m.Tire.DeclQuants =		MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.Tire.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_Tire, ThisModelKind_STI, &m);
}
