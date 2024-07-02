#INFOFILE1.1 (UTF-8) - Do not remove this line!
FileIdent = CarMaker-TestSeries 12
FileCreator = CarMaker 12.0.1
Description:
	TestAutomation to test the influence of the weight of the driver (Named Value Driver_wt) on the acceleration.
LastChange = 2023-09-28 12:12:36 rop
StartTime = 2022-10-20 14:38:08
EndTime = 2022-10-20 14:38:29
ReportTemplate =
Step.0 = Settings
Step.0.Name = Global Settings
Step.1 = TestRun
Step.1.Name = Tutorial/Acceleration_TestManager
Step.1.Param.0 = Driver_Wt NValue
Step.1.Char.0.Name = Throttle change to 80%
Step.1.Char.0.Description:
Step.1.Char.0.Identifier = Time
Step.1.Char.0.Unit =
Step.1.Char.0.Param.0 = RTexpr {DM.Long.StepTime >=4 ? DM.Gas= 0.8}
Step.1.Char.0.Param.1 = RTexpr {first() ? Qu::AccTime=0}
Step.1.Char.0.Param.2 = RTexpr {first (Car.Road.sRoad>=75) ? AccTime = DM.Long.StepTime}
Step.1.Crit.0.Name = Acceleration Time
Step.1.Crit.0.Description:
Step.1.Crit.0.Good = [get AccTime] <= 6.55
Step.1.Crit.0.Warn = [get AccTime] > 6.55 && [get AccTime] <= 6.65
Step.1.Crit.0.Bad = [get AccTime] > 6.65
Step.1.Diag.0.Name = Longitudinal Acceleration Vs Time
Step.1.Diag.0.Type = Line Diagram
Step.1.Diag.0.Mode = Quantity vs Time
Step.1.Diag.0.VarOpt = 1
Step.1.Diag.0.Grid = None
Step.1.Diag.0.NAxes = 1
Step.1.Diag.0.RefFile =
Step.1.Diag.0.Pic = SimOutput/OMEN/20221020/diagrams/Tutorial_Acceleration_TestManager_143828_2-8.png
Step.1.Diag.0.Param1.0 = Auto {} {} Time
Step.1.Diag.0.Param1.1 = Auto {} {} LongAcc
Step.1.Diag.0.Param2.0 = Time {} {} {}
Step.1.Diag.0.Param2.1 = Car.ax {} {} {}
Step.1.Diag.0.Param3.0 = {} {} {} {}
Step.1.Diag.0.Param3.1 = {} {} {} {}
Step.1.Diag.0.Param3.2 = {} {} {} {}
Step.1.Var.0.Name = Variation 0
Step.1.Var.0.Param = 70
Step.1.Var.0.Result = good
Step.1.Var.0.ResDate = 1666273095
Step.1.Var.0.ResFiles = SimOutput/OMEN/20221020/Tutorial_Acceleration_TestManager_143813.erg
Step.1.Var.0.ManLst = 0:long0 0:lat0
Step.1.Var.0.Char.0.Ref = Time
Step.1.Var.0.Char.0.Value = 6.823000000000613
Step.1.Var.0.Crit.0.Ref = Acceleration Time
Step.1.Var.0.Crit.0.Result = good
Step.1.Var.0.Diag.0.Pic =
Step.1.Var.1.Name = Variation 1
Step.1.Var.1.Param = 85
Step.1.Var.1.Result = warn
Step.1.Var.1.ResDate = 1666273102
Step.1.Var.1.ResFiles = SimOutput/OMEN/20221020/Tutorial_Acceleration_TestManager_143821.erg
Step.1.Var.1.ManLst = 0:long0 0:lat0
Step.1.Var.1.Char.0.Ref = Time
Step.1.Var.1.Char.0.Value = 6.916000000000644
Step.1.Var.1.Crit.0.Ref = Acceleration Time
Step.1.Var.1.Crit.0.Result = warn
Step.1.Var.1.Diag.0.Pic =
Step.1.Var.2.Name = Variation 2
Step.1.Var.2.Param = 100
Step.1.Var.2.Result = bad
Step.1.Var.2.ResDate = 1666273109
Step.1.Var.2.ResFiles = SimOutput/OMEN/20221020/Tutorial_Acceleration_TestManager_143828.erg
Step.1.Var.2.ManLst = 0:long0 0:lat0
Step.1.Var.2.Char.0.Ref = Time
Step.1.Var.2.Char.0.Value = 7.007000000000675
Step.1.Var.2.Crit.0.Ref = Acceleration Time
Step.1.Var.2.Crit.0.Result = bad
Step.1.Var.2.Diag.0.Pic =
