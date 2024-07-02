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
** Simple external suspension Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    int Susp_KnC_Register_MyModel (void);
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"
#include "MathUtils.h"

static const char ThisModelClass[]  	= "SuspKnC";
static const char ThisModelKind[]   	= "MyModel_F";
static const char ThisModelKindLR[] 	= "MyModel_FLR";
static const char ThisModelKindFrcCpl[]	= "MyModel_FrcCpl_R";
static const int  ThisVersionId     	= 1;

/* Model parameters (static) */
typedef struct tMyModel {
    double	qComp2tz;
    double	qSteer2rz;
    /* From here on only relevant for ForceCoupling */
    struct 	WC {
	double mass;		/* Mass Wheel Carrier */
	double MountSpring[3];	/* Mount point of Spring on Wheel Carrier */
	double MountDamper[3];
	double t[3];		/* DOFs of WC in Fr1 */
	double v[3];
	double a[3];
	double t_0[3];		/* DOFs of WC in Fr0 */
	double v_0[3];
	double a_0[3];
	double Tr2Fr0[3][3];
    }WC;
    double TireNomRadius;
    struct	Vehicle {
	double 	Joint2WC[3];	/* Initial position of Wheel Carrier */
	double	MountSpring[3]; /* Mount point of Spring on Chassis */
	double	MountDamper[3];
    }Vehicle;
} tMyModel;

static void
MyModel_Kinematics(void *MP, tSuspIF *IFMain, tSuspIF *IFOpp, double dt)
{
    /* In: 	Position of WC
     * Out: 	Spring/Damper TransAx
     * 		IF->Kin
     * 		IF->dKindq
     */
    struct tMyModel *mp = (struct tMyModel *)MP;
    double vecSpring[3], vecDamper[3], vec[3];

    /* Wheel Carrier */
    IFMain->Kin[ixtz] = mp->WC.t[2] - mp->Vehicle.Joint2WC[2];

    /* Spring */
    VEC_Add(vec, mp->WC.t, mp->WC.MountSpring);
    VEC_Sub(vecSpring, mp->Vehicle.MountSpring, vec);
    IFMain->Kin[ixSpring] = VEC_Norm(vecSpring);
    if (vecSpring[2]<0.0)
	IFMain->Kin[ixSpring] *= -1.0;

    /* Damper/Buffer/Stabi */
    VEC_Add(vec, mp->WC.t, mp->WC.MountDamper);
    VEC_Sub(vecDamper, mp->Vehicle.MountDamper, vec);
    IFMain->Kin[ixDamp]   = VEC_Norm(vecDamper);
    if (vecDamper[2]<0.0)
	IFMain->Kin[ixDamp] *= -1.0;
    IFMain->Kin[ixBuf]    = IFMain->Kin[ixDamp];
    IFMain->Kin[ixStabi]  = IFMain->qComp;

    /* Direction of action and point of attack of Force Elements */
    VEC_Normalize(IFMain->FrcCpl.Spring.TransAx, vecSpring);
    VEC_Normalize(IFMain->FrcCpl.Damper.TransAx, vecDamper);
    VEC_Assign(IFMain->FrcCpl.Spring.xtof, mp->Vehicle.MountSpring); /* const. */
    VEC_Assign(IFMain->FrcCpl.Damper.xtof, mp->Vehicle.MountDamper);

    /* Kinematics Gradients */
    IFMain->dqComp[ixtz] = 1.0;
    IFMain->dqComp[ixSpring] = -1.0/IFMain->Kin[ixSpring]*vecSpring[2];
    IFMain->dqComp[ixDamp] = -1.0/IFMain->Kin[ixDamp]*vecDamper[2];
    IFMain->dqComp[ixBuf] = IFMain->dqComp[ixDamp];
    IFMain->dqComp[ixStabi] = 1.0;
}

static void
MyModel_qCompCalc(void *MP, tSuspIF *IFMain, tSuspIF *IFOpp)
{
    /* In: 	Position AF
     * 		Position Tire Contact Point on Road
     * Out: 	IFMain->qComp and qpcomp
     */
    struct tMyModel *mp = (struct tMyModel *)MP;
    double vec[3], vec2[3];
    double P_0[3];

    VEC_Sub(P_0, IFMain->FrcCpl.tTireIF->P_0, IFMain->FrcCpl.AF.t_0);	/* AF to Tire_P_0 in Fr0 */
    VEC_MatTVec(vec2, IFMain->FrcCpl.AF.Tr2Fr0, P_0);			/* AF to Tire_P_0 in Fr1 */
    IFMain->qComp = mp->TireNomRadius - (-vec2[2] + mp->Vehicle.Joint2WC[2]);

    VEC_Cross(vec, IFMain->FrcCpl.AF.omega_0, P_0);
    VEC_AddS(vec, IFMain->FrcCpl.AF.v_0);
    VEC_MatTVec(vec2, IFMain->FrcCpl.AF.Tr2Fr0, vec);
    IFMain->qpComp = -vec[2];
}

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


static void
MyModel_Delete (void *MP)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Park the dict entries for dynamically allocated variables before deleting */
    MyModel_DeclQuants_dyn (mp, 1);
    free (mp);
}


static void *
MyModel_New (
    struct tInfos  	*Inf,
    struct tSuspCfgIF	*SuspCfgIF,
    const char     	*KindKey,
    const char     	*Pre)
{
    struct tMyModel *mp = NULL;
    char MsgPre[64];
    const char *ModelKind;
    int VersionId = 0;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_SuspKnC, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    mp = (struct tMyModel*)calloc(1, sizeof(*mp));

    /* SuspCfgIF */
    SuspCfgIF->Use_qSteer = 1;

    /* Kinematics */
    mp->qComp2tz  = iGetDblOpt (Inf, "MySusp_F.qComp2tz",  1.0);
    mp->qSteer2rz = iGetDblOpt (Inf, "MySusp_F.qSteer2rz", 5.0);

    return mp;
}


static void *
MyModel_New_LR (
    struct tInfos  	*Inf,
    struct tSuspCfgIF	*SuspCfgIF,
    const char     	*KindKey,
    const char     	*Pre)
{
    struct tMyModel *mp = NULL;
    char MsgPre[64];
    const char *ModelKind;
    int VersionId = 0;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_SuspKnC, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    mp = (struct tMyModel*)calloc(1, sizeof(*mp));

    /* SuspCfgIF */
    SuspCfgIF->Use_qSteer  = 1;
    SuspCfgIF->DoubleSided = 1;

    /* Kinematics */
    mp->qComp2tz  = iGetDblOpt (Inf, "MySusp_FLR.qComp2tz",  1.0);
    mp->qSteer2rz = iGetDblOpt (Inf, "MySusp_FLR.qSteer2rz", 5.0);

    return mp;
}

static void *
MyModel_New_FrcCpl (
	struct tInfos  	*Inf,
	struct tSuspCfgIF	*SuspCfgIF,
	const char     	*KindKey,
	const char     	*Pre)
{
    struct tMyModel *mp = NULL;
    struct tInfos *TireInf = SuspCfgIF->TireInf;
    char MsgPre[64];
    const char *ModelKind;
    int VersionId = 0;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_SuspKnC, KindKey,
					 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    mp = (struct tMyModel*)calloc(1, sizeof(*mp));

    /* SuspCfgIF */
    SuspCfgIF->DoubleSided = 0;
    SuspCfgIF->FrcFrame = ForceFrame_Fr1;
    SuspCfgIF->Use_qSteer = 0;
    SuspCfgIF->ForceCoupling = 1;

    /* Kinematics */
    mp->qComp2tz  = iGetDblOpt (Inf, "MySusp_F.qComp2tz",  1.0);
    mp->qSteer2rz = iGetDblOpt (Inf, "MySusp_F.qSteer2rz", 0.0); /* not used */

    /* Vehicle */
    VEC_Assign(mp->Vehicle.Joint2WC, SuspCfgIF->WheelBdy->pos);

    VEC_Assign(mp->Vehicle.MountSpring, mp->Vehicle.Joint2WC);
    mp->Vehicle.MountSpring[1] *= 0.95; /* inclined inward */
    mp->Vehicle.MountSpring[2] += 0.24; /* above WC */

    VEC_Assign(mp->Vehicle.MountDamper, mp->Vehicle.MountSpring);

    /* Wheel Carrier */
    mp->WC.mass = 25.0;

    VEC_Assign(mp->WC.MountSpring, Null3x1);
    VEC_Assign(mp->WC.MountDamper, Null3x1);

    VEC_Assign(mp->WC.t, mp->Vehicle.Joint2WC); // Initial position
    VEC_Assign(mp->WC.v, Null3x1); // Initial position

    /* Tire */
    mp->TireNomRadius = iGetDbl    (TireInf, "NomRadius");

    return mp;
}

static int
MyModel_Calc(void *MP, tSuspIF *IFMain, tSuspIF *IFOpp, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    int i;

    /* Kinematics on Main Side */
    IFMain->dqComp[ixtz] = 	 mp->qComp2tz;
    IFMain->dqComp[ixSpring] = 	-mp->qComp2tz;
    IFMain->dqComp[ixDamp] = 	-mp->qComp2tz;
    IFMain->dqComp[ixBuf] = 	-mp->qComp2tz;
    IFMain->dqComp[ixStabi] = 	 mp->qComp2tz;
    IFMain->dqSteer[ixrz] =	 mp->qSteer2rz;

    for (i=0; i<ixKinMax; i++) {
	IFMain->Kin[i] = IFMain->dqComp[i]  * IFMain->qComp
		       + IFMain->dqSteer[i] * IFMain->qSteer;
    }

    return 0;
}


static int
MyModel_Calc_LR (void *MP, tSuspIF *IFMain, tSuspIF *IFOpp, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    int i;

    /* Kinematics on Left Side */
    IFMain->dqComp[ixtz] = 	 mp->qComp2tz;
    IFMain->dqComp[ixSpring] = 	-mp->qComp2tz;
    IFMain->dqComp[ixDamp] = 	-mp->qComp2tz;
    IFMain->dqComp[ixBuf] = 	-mp->qComp2tz;
    IFMain->dqComp[ixStabi] = 	 mp->qComp2tz;
    IFMain->dqSteer[ixrz] =	 mp->qSteer2rz;

    for (i=0; i<ixKinMax; i++) {
	IFMain->Kin[i] = IFMain->dqComp[i]  * IFMain->qComp
		       + IFMain->dqSteer[i] * IFMain->qSteer;
    }

    /* Kinematics on Right Side */
    IFOpp->dqComp[ixtz] = 	 mp->qComp2tz;
    IFOpp->dqComp[ixSpring] = 	-mp->qComp2tz;
    IFOpp->dqComp[ixDamp] = 	-mp->qComp2tz;
    IFOpp->dqComp[ixBuf] = 	-mp->qComp2tz;
    IFOpp->dqComp[ixStabi] = 	 mp->qComp2tz;
    IFOpp->dqSteer[ixrz] =	 mp->qSteer2rz;

    for (i=0; i<ixKinMax; i++) {
	IFOpp->Kin[i] = IFOpp->dqComp[i]  * IFOpp->qComp
		      + IFOpp->dqSteer[i] * IFOpp->qSteer;
    }

    return 0;
}


static int
MyModel_Calc_FrcCpl(void *MP, tSuspIF *IFMain, tSuspIF *IFOpp, double dt)
{
    /* Model with 1 DOF and 5 Kinematic Constraints */
    struct tMyModel *mp = (struct tMyModel *)MP;
    double Frc[3], Trq[3];
    double t_0[3], v_0[3];
    double TrqCorr[3];

    /* Wheel compression */
    if (IFMain->FrcCpl.qCalc) {
	MyModel_qCompCalc(mp, IFMain, IFOpp);
    }

    /* Preparation phase/model check */
    if (dt == 0 || AppStartInfo.ModelCheck) {
	mp->WC.t[2] = mp->Vehicle.Joint2WC[2] + IFMain->qComp;
	mp->WC.v[2] = IFMain->qpComp;
    } /* else: mp->WC.t is inner state of model */

    /** Kinematics */
    MyModel_Kinematics(mp, IFMain, IFOpp, dt);

    /* Substitution models for Spring and Tire Force */
    if (dt == 0 || AppStartInfo.ModelCheck) {
	SuspMod_FrcSpring_Calc (IFMain->CfgIF->SuspNo,   IFMain->Kin[ixSpring], 0,
				&IFMain->FrcSpring);
	VEC_Mul(IFMain->Frc2WC, IFMain->FrcCpl.Spring.TransAx, IFMain->FrcSpring);
	IFMain->FrcDamp = 0.0;	/* static calculation */
	IFMain->FrcBuf = 0.0;	/* not considered in modelcheck and preparation phase */

    }

    /* Force Elements */
    VEC_Mul(Frc, IFMain->FrcCpl.Spring.TransAx, -IFMain->FrcSpring);
    VEC_AddMul(Frc, Frc, IFMain->FrcCpl.Damper.TransAx, -IFMain->FrcDamp);
    VEC_AddMul(Frc, Frc, IFMain->FrcCpl.Damper.TransAx, -IFMain->FrcBuf);
    /*  Tire Force/Torque */
    VEC_AddS(Frc, IFMain->Frc2WC);
    VEC_Assign(Trq, IFMain->Trq2WC);

    /** Dynamics */
    if (dt!=0.0) {
	mp->WC.a[2] = 1 / mp->WC.mass * Frc[2];
	/* Integration Euler Explicit */
	mp->WC.v[2] += mp->WC.a[2] * dt;
	mp->WC.t[2] += mp->WC.v[2] * dt;
    }

    /* Kinematics Update after Integration */
    MyModel_Kinematics(mp, IFMain, IFOpp, dt);

    /** Additional Model Output */
    /* Wheel Carrier global position */
    VEC_MatVec(t_0, IFMain->FrcCpl.AF.Tr2Fr0, mp->WC.t);
    VEC_MatVec(v_0, IFMain->FrcCpl.AF.Tr2Fr0, mp->WC.v);

    VEC_Add(IFMain->FrcCpl.WC.t_0, IFMain->FrcCpl.AF.t_0, t_0);
    VEC_Add(IFMain->FrcCpl.WC.v_0, IFMain->FrcCpl.AF.v_0, v_0);

    VEC_MatAssign(IFMain->FrcCpl.WC.Tr2Fr0, IFMain->FrcCpl.AF.Tr2Fr0);
    VEC_Assign(IFMain->FrcCpl.WC.omega_0, IFMain->FrcCpl.AF.omega_0);
    VEC_Assign(IFMain->FrcCpl.WC.alpha_0, IFMain->FrcCpl.AF.alpha_0);

    /* Forces to AF (inertia of Wheel Carrier is ignored) */
    Frc[2] = 0.0; /* no Forces in vertical direction (=DOF) */
    VEC_Assign(IFMain->FrcCpl.AF.Frc_1, Frc);
    VEC_Cross(TrqCorr, mp->WC.t, Frc); /* Correction Moment to Fr1*/
    VEC_Add(IFMain->FrcCpl.AF.Trq_1, TrqCorr, Trq);

    return 0;
}


int
Susp_KnC_Register_MyModel (void)
{
    tModelClassDescr m;

    memset (&m, 0, sizeof(m));
    m.SuspKnC.VersionId =	ThisVersionId;
    m.SuspKnC.New =		MyModel_New;
    m.SuspKnC.Calc =		MyModel_Calc;
    m.SuspKnC.Delete =		MyModel_Delete;
    m.SuspKnC.DeclQuants =	MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.SuspKnC.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_SuspKnC, ThisModelKind, &m);
}

int
Susp_KnC_Register_MyModel_LR (void)
{
    tModelClassDescr m;

    memset (&m, 0, sizeof(m));
    m.SuspKnC.VersionId =	ThisVersionId;
    m.SuspKnC.New =		MyModel_New_LR;
    m.SuspKnC.Calc =		MyModel_Calc_LR;
    m.SuspKnC.Delete =		MyModel_Delete;
    m.SuspKnC.DeclQuants =	MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.SuspKnC.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_SuspKnC, ThisModelKindLR, &m);
}

int
Susp_KnC_Register_MyModel_FrcCpl (void)
{
    tModelClassDescr m;

    memset (&m, 0, sizeof(m));
    m.SuspKnC.VersionId =	ThisVersionId;
    m.SuspKnC.New =		MyModel_New_FrcCpl;
    m.SuspKnC.Calc =		MyModel_Calc_FrcCpl;
    m.SuspKnC.Delete =		MyModel_Delete;
    m.SuspKnC.DeclQuants =	MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.SuspKnC.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_SuspKnC, ThisModelKindFrcCpl, &m);
}
