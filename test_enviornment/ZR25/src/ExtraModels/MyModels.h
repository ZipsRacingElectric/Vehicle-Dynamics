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

#ifndef _MYMODELS_H__
#define _MYMODELS_H__

#ifdef __cplusplus
extern "C" {
#endif

int VC_Register_AccelCtrl			(void);
int Aero_Register_MyModel			(void);
int Battery_Register_MyModel			(void);
int BatteryCU_Register_MyModel			(void);
int Brake_Register_MyModel			(void);
int Brake_Register_MyModelHydESPWrap		(void);
int Vehicle_Register_MyModel			(void);
int Clutch_Register_MyModel			(void);
int DriveLine_Register_MyModel			(void);
int DriveLineXWD_Register_MyModel		(void);
int Engine_Register_MyModel			(void);
int EngineCU_Register_MyModel			(void);
int Environment_Register_MyModel		(void);
int GearBox_Register_MyModel			(void);
int HydBrakeCU_Register_MyModel			(void);
int Motor_Register_MyModel			(void);
int MotorCU_Register_MyModel			(void);
int PowerSupply_Register_MyModel		(void);
int PowerTrain_Register_MyModel			(void);
int PowerTrainXWD_Register_MyModel		(void);
int PTControl_Register_MyModel			(void);
int PTControlOSM_Register_MyModel		(void);
int PTGenCoupling_Register_MyModel		(void);
int Steering_Register_MyModel			(void);
int Susp_Buffer_Register_MyModel		(void);
int Susp_Damper_Register_MyModel		(void);
int Susp_Spring_Register_MyModel		(void);
int Susp_ParasiticFriction_Register_MyModel	(void);
int Susp_ParasiticStiffness_Register_MyModel	(void);
int Susp_Stabi_Register_MyModel			(void);
int SuspEF_Buffer_Register_MyModel		(void);
int SuspEF_Damper_Register_MyModel		(void);
int SuspEF_Spring_Register_MyModel		(void);
int SuspEF_Stabi_Register_MyModel		(void);
int SuspExtFrcs_Register_MyModel		(void);
int Susp_KnC_Register_MyModel			(void);
int Susp_KnC_Register_MyModel_LR		(void);
int Susp_KnC_Register_MyModel_FrcCpl		(void);
int Susp_KnC_Register_MyModel_AuxFMU		(void);
int Tire_Register_MyModelCPI			(void);
int Tire_Register_MyModelSTI			(void);
int TireCPMod_Register_MyModel			(void);
int TireCPMod_Register_MyModelFourPoster	(void);
int Tire_Register_MyModelRTTireWrap		(void);
int TransmCU_Register_MyModel			(void);
int VehicleControl_Register_MyModel		(void);
int VhclOperator_Register_MyModel		(void);
int UserDriver_Register_MyModel			(void);

#ifdef __cplusplus
}
#endif

#endif	/* #ifndef _MYMODELS_H__ */
