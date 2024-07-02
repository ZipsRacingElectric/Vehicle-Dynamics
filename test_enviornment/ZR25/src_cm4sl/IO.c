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
** Connection to I/O hardware of the CarMaker/HIL test stand
**
** Connected test rig: ???
**
******************************************************************************
**
** Functions
** ---------
**
** - iGetCal ()
** - CalIn ()
** - CalInF ()
** - CalOut ()
** - CalOutF ()
** - LimitInt ()
** - IO_Init_First ()
** - IO_Init_Finalize ()
** - IO_Init ()
** - IO_Param_Get ()
** - IO_BeginCycle ()
** - IO_In ()
** - IO_Out ()
** - IO_Cleanup ()
**
******************************************************************************
*/

#include <Global.h>

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

#include <CarMaker.h>

#if defined(XENO)
#include <mio.h>
#endif /* defined(XENO) */
#include <ioconf.h>
#if defined(CM_HIL)
#include <FailSafeTester.h>
#endif /* defined(CM_HIL) */

#include "IOVec.h"


/*** I/O vector */
tIOVec IO;


/*** I/O configuration */

/* int IO_None; DON'T - Variable is predefined by CarMaker! */
int IO_CAN_IF;
int IO_FlexRay;

static struct tIOConfig IOConfiguration[] = {
/* This table should contain one line for each IO_xyz-flag in IOVec.h */
/*  { <Flagvar>,	<Name for -io>,	<Description for -help> },     */
    { &IO_None,		"none",		"No I/O" }, /* Always keep this first line! */
    { &IO_CAN_IF,	"can",		"CAN communication" },
    { &IO_FlexRay,	"flexray",	"FlexRay communication" },
    { NULL, NULL, NULL } /* End of table */
};



/**** Additional useful functions *********************************************/


/*
** iGetCal()
**
** Read calibration parameters.
*/

void
iGetCal (tInfos *Inf, const char *key, tCal *cal, int optional)
{
    cal->Min       =  1e37f;
    cal->Max       = -1e37f;
    cal->LimitLow  =  1e37f;
    cal->LimitHigh = -1e37f;
    cal->Factor    =  1.0f;
    cal->Offset    =  0.0f;
    cal->Rezip     =  0;

    const char *item = iGetStrOpt(Inf, key, NULL);
    if (item != NULL) {
	int n = sscanf(item, "%g %g %g %g %d",
		       &cal->LimitLow, &cal->LimitHigh,
		       &cal->Factor,   &cal->Offset, &cal->Rezip);
	if (n != 5)
	    LogErrF(EC_Init, "Invalid calibration parameter entry '%s'", key);
    } else {
	if (!optional) {
	    LogErrF(EC_Init, "Missing calibration parameter entry '%s'", key);
	    return;
	}
	cal->LimitLow  = -1e37f;
	cal->LimitHigh =  1e37f;
    }

    cal->Min = cal->LimitHigh;
    cal->Max = cal->LimitLow;
}


/*
** CalInF() / CalIn()
**
** Analog input -> calibration infos -> physical quantity
** Converts an I/O value (e.g. the voltage from an analog input module) to
** the corresponding physical value, delimited by LimitLow and LimitHigh.
*/

float
CalInF (tCal *cal, float Value)
{
    float Result = (Value - cal->Offset) * cal->Factor;

    if (cal->Rezip)
	Result = 1.0f / Result;

    if      (Result < cal->Min)  cal->Min = Result;
    else if (Result > cal->Max)  cal->Max = Result;

    if      (Result < cal->LimitLow)   Result = cal->LimitLow;
    else if (Result > cal->LimitHigh)  Result = cal->LimitHigh;

    return Result;
}

float
CalIn (tCal *cal, int Value)
{
    return CalInF(cal, (float) Value);
}


/*
** CalOutF() / CalOut()
**
** Physical quantity -> calibration infos -> analog output
** The physical value is delimited by LimitLow and LimitHigh and then converted
** to the corresponding I/O value (e.g. voltage for an analog output module).
*/

float
CalOutF (tCal *cal, float Value)
{
    if      (Value < cal->Min) cal->Min = Value;
    else if (Value > cal->Max) cal->Max = Value;

    if      (Value < cal->LimitLow)  Value = cal->LimitLow;
    else if (Value > cal->LimitHigh) Value = cal->LimitHigh;

    if (cal->Rezip) {
	return 1.0f / (Value*cal->Factor) + cal->Offset;
    } else {
	return Value/cal->Factor + cal->Offset;
    }
}

int
CalOut (tCal *cal, float Value)
{
    return (int)CalOutF(cal, Value);
}


int
LimitInt (float fValue, int Min, int Max)
{
    int Value = (int)fValue;
    if      (Value < Min) return Min;
    else if (Value > Max) return Max;
    return   Value;
}



/*****************************************************************************/


/*
** IO_Init_First ()
**
** First, low level initialization of the IO module
**
** Call:
** - one times at start of program
** - no realtime conditions
*/

int
IO_Init_First (void)
{
    memset(&IO,  0, sizeof(IO));

    IO_SetConfigurations(IOConfiguration);

    return 0;
}



/*
** IO_Init ()
**
** initialization
** - i/o hardware
** - add variables to data dictionary
**
** call:
** - single call at program start
*/

int
IO_Init (void)
{
    Log("I/O Configuration: %s\n", IO_ListNames(NULL, 1));

#if defined(XENO)
    if (IOConf_Init() < 0)
            return -1;
#endif /* defined(XENO) */

    /* hardware configuration "none" */
    if (IO_None)
	return 0;
#if defined(CM_HIL)
    int nErrors = Log_nError;

    /*** MIO initialization */
    if (MIO_Init(NULL) < 0) {
	LogErrF(EC_General, "MIO initialization failed. I/O disabled (1)");
	IO_SelectNone();
	return -1;
    }
    MIO_SetAppState(0.0, MIO_SimState_AppInit);
    // MIO_ModuleShow();

    /* check for MIO errors */
    if (nErrors != Log_nError) {
	LogErrF(EC_General, "MIO initization failed. I/O disabled (2)");
	IO_SelectNone();
	return -1;
    }
#endif /* defined(CM_HIL) */

#if defined(XENO)
    /*** FailSafeTester */
    FST_ConfigureCAN();
#endif /* defined(XENO) */

    return 0;
}


/*
** IO_Init_Finalize ()
**
** last (deferred) I/O initialization step
**
** call:
** - single call at program start in CarMaker_FinishStartup()
*/

int
IO_Init_Finalize (void)
{

    return 0;
}


/*
** IO_Param_Get ()
**
** Get i/o configuration parameters
** - calibration
** - constant values
** - ids
*/

int
IO_Param_Get (tInfos *inf)
{
    unsigned nError = GetInfoErrorCount ();
#if defined(CM_HIL)
    /* ignition off */
    SetKl15 (0);
#if defined(XENO)
    IOConf_Param_Get();
#endif /* defined(XENO) */
#endif /* defined(CM_HIL) */

    if (IO_None)
    	return 0;

    return nError != GetInfoErrorCount() ? -1 : 0;
}


void
IO_BeginCycle (void)
{
#if defined(CM_HIL)
    MIO_SetAppState(TimeGlobal, (tMIO_SimState)SimCore_State2MIO_SimState(SimCore.State));
#endif /* defined(CM_HIL) */
}


/*
** IO_In ()
**
** reading signals from hardware / ECU
**
** CycleNo: simulation cycle counter, incremented every loop/millisecond
**
** call:
** - in the main loop
** - first function call in main loop, after waiting for next loop
** - just before User_In()
** - pay attention to realtime condition
*/

void
IO_In (unsigned CycleNo)
{
#if defined(CM_HIL)
    CAN_Msg Msg;

    IO.DeltaT = SimCore.DeltaT;
    IO.T      = TimeGlobal;
#if defined(XENO)
    IOConf_In(CycleNo);
#endif /* defined(XENO) */
    if (IO_None)
	return;

    /*** FailSafeTester messages */
    if (FST_IsActive()) {
#if defined(XENO)
	while (MIO_M51_Recv(FST_CAN_Slot, FST_CAN_Ch, &Msg) == 0)
	    FST_MsgIn (CycleNo, &Msg);
#endif /* defined(XENO) */
    }
#endif /* defined(CM_HIL) */

}



/*
** IO_Out ()
**
** writing signals to hardware / ECU
**
** CycleNo: simulation cycle counter, incremented every loop/millisecond
**
** call:
** - in the main loop
** - last function call in main loop
** - just after User_Out()
** - pay attention to realtime condition
*/

void
IO_Out (unsigned CycleNo)
{
#if defined(XENO)
    IOConf_Out(CycleNo);
#endif /* defined(XENO) */
    if (IO_None)
	return;
#if defined(CM_HIL)
    /*** Messages to the FailSafeTester */
    FST_MsgOut(CycleNo);
#endif /* defined(CM_HIL) */

}



/*
** IO_Cleanup ()
**
** Uninits all MIO hardware:
** - puts M-Modules into reset state
** - frees unneeded memory
*/

void
IO_Cleanup (void)
{
    if (IO_None)
	goto EndReturn;

#if defined(XENO)
    IOConf_Cleanup();
#endif
#if defined(CM_HIL)
    MIO_SetAppState(TimeGlobal, MIO_SimState_AppExit);
    MIO_ResetModules();
    MIO_DeleteAll();
#endif /* defined(CM_HIL) */

  EndReturn:
    return;
}
