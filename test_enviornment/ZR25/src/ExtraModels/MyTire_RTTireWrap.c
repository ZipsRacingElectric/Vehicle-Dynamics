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
** A Wrapper around an existing, already registered RTTire model
**
** Add the declaration of the register function to one of your header files,
** for example to User.h and call it in User_Register()
**
**	Tire_Register_MyModelRTTireWrap();
**
******************************************************************************
*/

#include <Global.h>

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "CarMaker.h"
#include "Car/Vehicle_Car.h"
#include "MyModels.h"

static const char ThisModelClass[] = "Tire";
static const char ThisModelKind[]  = "RTTireWrap";
static const int  ThisVersionId    = 1;

struct tMyModel {
    /* Parameters */
    tTParamHeader Head;			/* Struct must begin with this header */

    struct {
	char		*path;		/* Path of RTTire file */
	time_t		ModTime;	/* Modification time of RTTire file */
    } TireFile;

    /* wrapped model: functions New/Calc/Delete/... */
#ifdef __cplusplus
    struct tModelClassDescr::tModelClassDescr_Tire Wrapped;
#else
    struct tModelClassDescr_Tire                   Wrapped;
#endif

    void		*Wrapped_param;		/* wrapped model: 	      *
						 * parameter handle           */
    double		muFac;			/* muRoad factor	      */
};


static void
MyModel_DeclQuants_dyn (struct tMyModel *mp, int park)
{
    char buf[64];
    static struct tMyModel MyModel_Dummy = {{0}};
    if (park)
	mp = &MyModel_Dummy;

    /* Define here dict entries for dynamically allocated variables. */

    /* Call wrapped tire model quantities */
    if (mp->Wrapped.DeclQuants != NULL)
	mp->Wrapped.DeclQuants(mp->Wrapped_param);

    sprintf (buf, "RTTireWrap.muFac_%s", Vehicle_TireNo_Str(mp->Head.TireNo));
    DDefDouble4 (NULL, buf, "", &mp->muFac, DVA_VC);
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


static int
MyModel_ParamsChanged (void *MP)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    if (mp->TireFile.path==NULL ||
	mp->TireFile.ModTime != PathModTime(mp->TireFile.path))
	return 1;

    return 0;
}


static void *
MyModel_New (
    struct tInfos	*inf,
    const char		*ParamDir,
    int			TireNo,
    int			Side,
    int			Twin)
{
    struct tMyModel *mp = NULL;
    tModelClassDescr *md;
    const char 	*key;
    char        MsgPre[64], *path, ModelKind[128] = "";
    int		VersionId = 0;
    char	*item;
    const char *FName;
    struct tInfos *Inf_RTTire;

    sprintf (MsgPre, "%s %s", ThisModelClass, ThisModelKind);

    FName = ((FName=InfoGetFilename(inf)) != NULL ? FName : "(unknown)");

    /*** check ModelKind, ModelVersion */
    item = iGetStrOpt(inf, "FileIdent", "");
    if (sscanf(item, "CarMaker-Tire-%s %d", ModelKind, &VersionId) < 1) {
	LogErrF (EC_Init, "%s: Wrong syntax in tire kind "
		 "information (FileIdent=%s, file '%s')",
		 MsgPre, item, FName);
	goto ErrorReturn;
    }

    if (strcmp(ModelKind, ThisModelKind) != 0) {
	LogErrF (EC_Init,
		 "%s: Tire with wrong Kind"
		 " (got '%s', expected '%s %d', file '%s')",
		 MsgPre, item, ThisModelKind, ThisVersionId, FName);
	goto ErrorReturn;
    }
    if (VersionId != ThisVersionId) {
	LogErrF (EC_Init, "%s: Tire model '%s' "
		 "with wrong version (got %d, expected %d, file '%s')",
		 MsgPre, ThisModelKind, VersionId, ThisVersionId, FName);
	goto ErrorReturn;
    }

    mp = (struct tMyModel*)calloc(1, sizeof(struct tMyModel));
    mp->muFac = 1.0;

    key = "RTTire_FName";
    path = iGetStr(inf, "RTTire_FName");
    if (path == NULL || strlen(path) == 0) {
	LogErrF(EC_Init, "%s: Missing RTTire file name '%s'", MsgPre, key);
	goto ErrorReturn;
    }
    path = PathJoin(2, ParamDir, path);
    mp->TireFile.path    = strdup(path);
    mp->TireFile.ModTime = PathModTime(mp->TireFile.path);

    Inf_RTTire = InfoNew();
    if (iRead2(NULL, Inf_RTTire, path, MsgPre) < 0)
	goto ErrorReturn;

    /* Get interface functions of wrapped model. */
    key = "RTTire";
    if ((md = Model_Lookup(ModelClass_Tire, key)) == NULL) {
	LogErrF(EC_Init, "%s: Missing tire model '%s'", MsgPre, key);
	goto ErrorReturn;
    }
    mp->Wrapped = md->Tire;

    /* Initialize wrapped model. */
    if (mp->Wrapped.New != NULL) {
	mp->Wrapped_param = mp->Wrapped.New(Inf_RTTire, ParamDir, TireNo, Side, Twin);
	if (mp->Wrapped_param == NULL)
	    goto ErrorReturn;
	memcpy (&mp->Head, mp->Wrapped_param, sizeof(tTParamHeader));
    }

    return mp;

   ErrorReturn:
	free(mp);
	return NULL;
}


static int
MyModel_Calc (void *MP, tTireIF *IF, double dt)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Modify input variables */
    IF->muRoad *= mp->muFac;

    /* Call wrapped RTTire brake model. */
    if (mp->Wrapped.Calc (mp->Wrapped_param, IF, dt) < 0)
	return -1;

    /* Modify output variables */

    return 0;
}


static void
MyModel_Delete (void *MP)
{
    struct tMyModel *mp = (struct tMyModel *)MP;

    /* Park the dict entries for dynamically allocated variables before deleting */
    MyModel_DeclQuants_dyn (mp, 1);

    if (mp->Wrapped_param != NULL && mp->Wrapped.Delete != NULL)
	mp->Wrapped.Delete (mp->Wrapped_param);

    if (mp->TireFile.path != NULL)
	free (mp->TireFile.path);

    free (mp);
}


int 
Tire_Register_MyModelRTTireWrap (void)
{
    tModelClassDescr m;

    memset(&m, 0, sizeof(m));
    m.Tire.VersionId =		ThisVersionId;
    m.Tire.is3DTire = 		0;	// set 1 if the STI interface is used
    m.Tire.New =		MyModel_New;
    m.Tire.Calc =		MyModel_Calc;
    m.Tire.Delete =		MyModel_Delete;
    m.Tire.DeclQuants =		MyModel_DeclQuants;
    m.Tire.ParamsChanged = 	MyModel_ParamsChanged;

    return Model_Register(ModelClass_Tire, ThisModelKind, &m);
}
