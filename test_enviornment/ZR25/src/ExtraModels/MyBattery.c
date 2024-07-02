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
** Simple Battery Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    Battery_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "PowerTrain.PowerSupply.Batt";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;


struct tMyModel {
    char	Ident[20];
    double	Capacity;
    double	R0;
    double      Volt_oc0;
    double      Temp;
    double	Pwr_max;

    double	SOC;
    double	Volt_oc;
    double	Volt0;
};


static void
MyModel_DeclQuants_dyn (struct tMyModel *mp, int park)
{
    tDDefault *df = DDefaultCreate("PT.Batt%s.", mp->Ident);

    static struct tMyModel MyModel_Dummy = {{0}};
    if (park)
	mp = &MyModel_Dummy;

    DDefDouble4 (df, "Volt_oc",	"V", 	&mp->Volt_oc,       DVA_None);

    DDefaultDelete(df);

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
    struct tMyModel *mp = (struct tMyModel *) MP;

    /* Park the dict entries for dynamically allocated variables before deleting */
    MyModel_DeclQuants_dyn (mp, 1);
    free (mp);
}


/* Model output parameters in the configuration struct CfgIF, which are required
   by CarMaker, are read in before the MyModel_New() function.
   - The parametrization of these parameters is supported by the GUI.
   - These output parameters can be used internally by the model in same way like
     the input parameters
*/
static void *
MyModel_New (struct tInfos *Inf, struct tPTBatteryCfgIF *CfgIF, const char *KindKey, const char *Ident)
{
    struct tMyModel *mp = NULL;
    char *key, buf[64], MsgPre[64], PreKey[64];
    const char *ModelKind;
    int VersionId = 0;

    sprintf(MsgPre, "%s %s", ThisModelClass, ThisModelKind);
    sprintf(PreKey, "%s%s.", ThisModelClass, Ident);

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_PTBattery, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel*) calloc(1, sizeof(*mp));
    strcpy (mp->Ident, Ident);

    /* get CfgIF parameters */
    if (Battery_GetCfgOutIF(Inf, CfgIF, ModelKind, Ident) != 0)
	goto ErrorReturn;

    /* CfgIF output -> Model */
    mp->Capacity = CfgIF->Capacity;
    mp->Volt_oc0 = CfgIF->Voltage;

    /* Battery temperature */
    mp->Temp = CfgIF->Temp_init; 

    /* R0 */
    mp->R0 = fabs(iGetDblOpt(Inf, key = strcat(strcpy(buf, PreKey), "R0"), 0.0012));
    if (mp->R0 <= 0.0) {
	LogErrF(EC_Init, "%s: Parameter '%s' must be positive and non zero", MsgPre, key);
	goto ErrorReturn;
    }

    /* Pwr_max */
    mp->Pwr_max = iGetDblOpt(Inf, strcat(strcpy(buf, PreKey), "Pwr_max"), 100.0) * 1e3;

    /* Model Quantities */
    MyModel_DeclQuants(mp);

    return mp;

    ErrorReturn:
	if (mp != NULL)
	    MyModel_Delete(mp);
	return NULL;
}


static int
MyModel_Calc (void *MP, struct tPTBatteryIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Amount of Charge */
    IF->AOC -= IF->Current * dt / 3600.0;
    IF->AOC  = M_BOUND(0.0, mp->Capacity, IF->AOC);

    /* State of Charge */
    mp->SOC = IF->AOC / mp->Capacity * 100.0;

    /* Voltage */
    mp->Volt_oc = mp->Volt_oc0;
    if (mp->SOC <= 1e-2)
	mp->Volt_oc = 0.0;
    mp->Volt0  = IF->Current * mp->R0;

    IF->Voltage = mp->Volt_oc - mp->Volt0;
    IF->Voltage = M_MAX(IF->Voltage, 0.0);

    /* Energy in kWh */
    IF->Energy = IF->AOC * IF->Voltage * 1e-3;

    /* Temperature */
    IF->Temp         = mp->Temp;
    IF->TempCool_out = IF->TempCool_in;

    /* Pwr_max */
    IF->Pwr_max = mp->Pwr_max;

    return 0;
}


int 
Battery_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.PTBattery.VersionId =	ThisVersionId;
    m.PTBattery.New =		MyModel_New;
    m.PTBattery.Calc =		MyModel_Calc;
    m.PTBattery.Delete =	MyModel_Delete;
    m.PTBattery.DeclQuants =	MyModel_DeclQuants;
    /* Should only be used if the model doesn't read params from extra files */
    m.PTBattery.ParamsChanged = ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_PTBattery, ThisModelKind, &m);
}
