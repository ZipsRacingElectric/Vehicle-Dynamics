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
** Simple gearbox Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    GearBox_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "PowerTrain.GearBox";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


struct tMyModel {
    int		nFGears;	        /* Number of forward gears */
    double 	iF[VEHICLE_NGEARS+1];	/* Ratio table of forward gears */
    int		nBGears;		/* Number of backward gears */
    double 	iB[VEHICLE_NGEARS+1];	/* Ratio table of backward gears */
    double	I_in, I_out;		/* Inertia in and out */
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
MyModel_New (tInfos *Inf, struct tPTGearBoxCfgIF *CfgIF, const char *KindKey, const char *Ident)
{
    struct tMyModel *mp = NULL;
    char 	MsgPre[64];
    int		i, VersionId = 0;
    const char *ModelKind;

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_PTGearBox, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*)calloc(1,sizeof(*mp));

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    /* get CfgIF parameters */
    if (GearBox_GetCfgOutIF (Inf, CfgIF, ModelKind, Ident) != 0)
	goto ErrorReturn;

    /* Inertia In/Out */
    mp->I_in  = iGetDblOpt (Inf, "PowerTrain.GearBox.I_in",  1e-4);
    mp->I_out = iGetDblOpt (Inf, "PowerTrain.GearBox.I_out", 0.05);

    /* CfgIF output -> Model */
    mp->nFGears = CfgIF->nFGears;
    for (i=0; i <= mp->nFGears; i++)
	mp->iF[i] = CfgIF->iFGear[i];

    mp->nBGears = CfgIF->nBGears;
    for (i=0; i <= mp->nBGears; i++)
	mp->iB[i] = CfgIF->iBGear[i];

    /* CfgIF output: verification if the parametrization by the GUI corresponds to the model */
    if (CfgIF->GBKind != GBKind_Manual) {
	LogErrF (EC_Init, "%s: model supports only manual gearbox", MsgPre);
	goto ErrorReturn;
    }
    if (CfgIF->ClKind != ClKind_Closed) {
	LogErrF (EC_Init, "%s: model supports only 'closed' clutch", MsgPre);
	goto ErrorReturn;
    }

    return mp;

    ErrorReturn:
	free (mp);
	return NULL;
}


static int
MyModel_PreSimSetup (void *MP, struct tPTGearBoxPreSimIF *PreSimIF,
		      	       struct tPTGearBoxIF       *IF)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    IF->GearNo   = PreSimIF->GearNo;
    IF->rotv_out = PreSimIF->rotv_out;
    IF->i        = (IF->GearNo>0) ? mp->iF[IF->GearNo] : mp->iB[-IF->GearNo];
    IF->rotv_in  = IF->i * IF->rotv_out;

    return 0;
}


static int
MyModel_Calc (void *MP, struct tPTGearBoxIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    double Inert_in = IF->Inert_in + mp->I_in;

    IF->GearNo  =  IF->GearNoTrg;
    IF->i       = (IF->GearNo>0) ? mp->iF[IF->GearNo] : mp->iB[-IF->GearNo];
    IF->Trq_out = IF->i * IF->Trq_in;

    if (IF->i != 0.0) {
	IF->rotv_in       = IF->i * IF->rotv_out;
	IF->Inert_out     = Inert_in * IF->i * IF->i + mp->I_out;
	IF->Trq_SuppInert = 0.0;

    } else {
	/* Free rotation -> DOF */
	double rota       = IF->Trq_in / Inert_in;
	IF->rotv_in      += rota * dt;
	IF->Inert_out     = mp->I_out;
	IF->Trq_SuppInert = rota * Inert_in;
    }
    IF->rot_in += IF->rotv_in * dt;

    IF->i_TrqIn2Out = IF->i;

    return 0;
}


static int
MyModel_ModelCheck (void *MP, tInfos *Inf)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    int	 i, n;
    FILE *fp = NULL;

    if ((fp = ModelCheck_GetDesignFile()) == NULL)
	return -1;

    fprintf (fp, "### GearBox.Kind = %s\n", ThisModelKind);

    fprintf (fp, "GearBox.I_in =            %10.7f\n", mp->I_in);
    fprintf (fp, "GearBox.I_out =           %10.7f\n", mp->I_out);

    fprintf (fp, "Ratio forward  gears [%d]: ", mp->nFGears);
    for (i=1, n=mp->nFGears+1; i < n; i++) {
	fprintf (fp, "  %10.7f", mp->iF[i]);
    }
    fprintf (fp, "\n");

    fprintf (fp, "Ratio backward gears [%d]: ", mp->nBGears);
    for (i=1, n=mp->nBGears+1; i < n; i++) {
	fprintf (fp, "  %10.7f", mp->iB[i]);
    }
    fprintf (fp, "\n");
    fprintf (fp, "\n");

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
GearBox_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.PTGearBox.VersionId =	ThisVersionId;
    m.PTGearBox.New =		MyModel_New;
    m.PTGearBox.Calc =		MyModel_Calc;
    m.PTGearBox.Delete =	MyModel_Delete;
    m.PTGearBox.DeclQuants =	MyModel_DeclQuants;
    m.PTGearBox.ModelCheck =	MyModel_ModelCheck;
    m.PTGearBox.PreSimSetup =	MyModel_PreSimSetup;
    /* Should only be used if the model doesn't read params from extra files */
    m.PTGearBox.ParamsChanged = ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_PTGearBox, ThisModelKind, &m);
}
