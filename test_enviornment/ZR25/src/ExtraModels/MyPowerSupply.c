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
** Simple power supply Model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**    PowerSupply_Register_MyModel ();
**
******************************************************************************
*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "PowerTrain.PowerSupply";
static const char ThisModelKind[]  = "MyModel";
static const int  ThisVersionId    = 1;

typedef struct tMyBattery {
    double	Capacity;
    double	R0;
    double	Volt_oc0;
    double	Pwr_max;

    double	Current;
    double	Voltage;
    double	Energy;
    double	SOC;
    double	AOC;
    double	SOC_init;
    double	Volt_oc;
    double	Volt0;

    double	Temp;
} tMyBattery;

struct tMyModel {
    struct tMyBattery LV;
};

/* Model output parameters in the configuration struct CfgIF, which are required
   by CarMaker, are read in before the MyModel_New() function.
   - The parametrization of these parameters is supported by the GUI.
   - These output parameters can be used internally by the model in same way like
     the input parameters
*/
static void *
MyModel_New (struct tInfos *Inf, struct tPTPowerSupplyCfgIF *CfgIF, const char *KindKey)
{
    struct tMyModel *mp   = NULL;
    struct tMyBattery *LV = NULL;
    char MsgPre[64], *key, buf[64], PreKey[64];
    const char *ModelKind;
    int VersionId = 0;
    double SOC_init;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    if ((ModelKind = SimCore_GetKindInfo(Inf, ModelClass_PTPowerSupply, KindKey,
	 				 0, ThisVersionId, &VersionId)) == NULL)
	return NULL;

    mp = (struct tMyModel *)calloc(1,sizeof(*mp));

    /* get CfgIF parameters */
    if (PowerSupply_GetCfgOutIF (Inf, CfgIF, ModelKind) != 0)
	goto ErrorReturn;

    LV = &mp->LV;

    /* Creating a new Battery - see also MyBattery.c */
    sprintf (PreKey, "%s.BattLV.", ThisModelClass);

    /* R0 */
    LV->R0 = fabs(iGetDblOpt (Inf, key=strcat(strcpy(buf, PreKey), "R0"), 0.0012));
    if (LV->R0==0.0) {
	LogErrF (EC_Init, "%s: Parameter '%s' must be positive and non zero", MsgPre, key);
	goto ErrorReturn;
    }

    /* Pwr_max */
    LV->Pwr_max = iGetDblOpt (Inf, strcat(strcpy(buf, PreKey), "Pwr_max"), 100.0) * 1e3;

    /* CfgIF output -> Model */
    LV->Capacity = CfgIF->BattLV.Capacity;
    LV->Volt_oc0 = CfgIF->BattLV.Voltage;
    LV->Temp	 = CfgIF->BattLV.Temp_init;

    /* Initialization of AOC */
    SOC_init = iGetDblOpt (Inf, strcat(strcpy(buf, PreKey), "SOC"), 70.0);
    SOC_init = M_BOUND (0.0, 100.0, SOC_init);
    LV->AOC  = SOC_init / 100.0 * LV->Capacity;

    return mp;

    ErrorReturn:
	free (mp);
	return NULL;
}


static int
MyModel_Calc (void *MP, struct tPTPowerSupplyIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;
    struct tMyBattery *LV = &mp->LV;

    if (IF->Voltage_LV > 0.0)
	IF->BattLV.Current = IF->Pwr_LV / IF->Voltage_LV;
    else
	IF->BattLV.Current = 0.0;

    LV->Current = IF->BattLV.Current;

    /* Battery Calculation - see also MyBattery.c */

    LV->AOC -= LV->Current * dt / 3600.0;
    LV->AOC  = M_BOUND (0.0, LV->Capacity, LV->AOC);

    LV->SOC = LV->AOC / LV->Capacity * 100.0;

    LV->Volt_oc = LV->Volt_oc0;
    if (LV->SOC <= 1e-2)
	LV->Volt_oc = 0.0;

    LV->Volt0  = LV->Current * LV->R0;

    LV->Voltage = LV->Volt_oc - LV->Volt0;
    LV->Voltage = M_MAX(LV->Voltage, 0.0);

    LV->Energy = LV->AOC * LV->Voltage * 1e-3;

    IF->BattLV.Pwr_max = LV->Pwr_max;

    /* Temperature */
    IF->BattLV.Temp = LV->Temp;
    IF->BattLV.TempCool_out = IF->BattLV.TempCool_in;

    IF->Voltage_LV     = LV->Voltage;
    IF->BattLV.AOC     = LV->AOC;
    IF->BattLV.Energy  = LV->Energy;

    return 0;
}


static void
MyModel_Delete (void *MP)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    free (mp);
}


int 
PowerSupply_Register_MyModel (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.PTPowerSupply.VersionId =	ThisVersionId;
    m.PTPowerSupply.New =		MyModel_New;
    m.PTPowerSupply.Calc =		MyModel_Calc;
    m.PTPowerSupply.Delete =		MyModel_Delete;
    m.PTPowerSupply.DeclQuants =	NULL;
    m.PTPowerSupply.ModelCheck =	NULL;
    /* Should only be used if the model doesn't read params from extra files */
    m.PTPowerSupply.ParamsChanged = 	ParamsChanged_IgnoreCheck;

    return Model_Register(ModelClass_PTPowerSupply, ThisModelKind, &m);
}
