#! /bin/sh
# next line restarts using wish \
exec CM -wish "$0" "$@"
#
# Instruments.tcl
#
# Copyright (c) IPG Automotive GmbH, Karlsruhe, Germany (www.ipg.de)

set Pgm(Dir)  [file dirname [info script]]

set Pgm(Name)           "Instruments"
set Pgm(Version)	"5.1"
set Pgm(VersionDate)    "2016-2-29"

array set Config {
    MaxSpeed 	220
    MaxRPM	7000
    MaxStW	550
    nGears	6
    Mode	CM
    Automatic   0
    AutoManual	0
    Units	""
    FuelEff	""
    ElecEff	""
    PTKind	"Generic"
}

# Subscription = list of variables to be subscribed
#   Tcl var			quantity        alternate quantities
set Subscription {
    {Qu(Distance) 		Vhcl.Distance	}
    {Qu(vCar) 			IO.vCar		Vhcl.v }
    {Qu(nEngine) 		IO.nEngine	Vhcl.Engine.rotv }
    {Qu(Gas) 			IO.Gas		DM.Gas		}
    {Qu(Clutch) 		IO.Clutch	DM.Clutch	}
    {Qu(Brake) 			DM.Brake	}
    {Qu(BrakeLever) 		DM.BrakeLever	}
    {Qu(GearNo) 		DM.GearNo	}
    {Qu(SelectorCtrl) 		DM.SelectorCtrl	}
    {Qu(PT.GearNo) 		PT.GearBox.GearNo DM.GearNo	}
    {Qu(StWhlAngle) 		Vhcl.Steer.Ang	}
    {Qu(SpeedLimit) 		DM.SpeedLimit	}
    {Qu(FuelConsump_Abs)	PT.Control.Consump.Fuel.Abs PT.ECU.Consump_Abs PT.Engine.Consump_Abs }
    {Qu(FuelConsump_Avg)	PT.Control.Consump.Fuel.Avg PT.ECU.Consump_Avg PT.Engine.Consump_Avg }
    {Qu(FuelConsump_Act)	PT.Control.Consump.Fuel.Act PT.ECU.Consump_Act PT.Engine.Consump_Act }
    {Qu(ElecConsump_Abs)	PT.Control.Consump.Elec.Abs  }
    {Qu(ElecConsump_Avg)	PT.Control.Consump.Elec.Avg  }
    {Qu(ElecConsump_Act)	PT.Control.Consump.Elec.Act  }
    {Qu(Temperature)		Ambient.Temperature }
    {Qu(LaABS) 			IO.LaABS	ESP.ABS.Active }
    {Qu(LaESP) 			IO.LaESP	IO.LaASRFDR	ESP.ESP.lamp }
    {Qu(_IsActive)		AccelCtrl.ACC.IsActive }
    {Qu(ACC_Targ_Dtct)		DASensor.RadarL.relvTgt.dtct }
    {Qu(ACC_DesiredSpd)		AccelCtrl.ACC.DesiredSpd }
    {Qu(EngineSwitch)		DM.SST		DM.EngineSwitch }
    {Qu(IndL)			VC.Lights.IndL }
    {Qu(IndR)			VC.Lights.IndR }
    {Qu(User.ESP_off)		ESP.User.ESP_off }
    {Qu(StratMode)		PT.Control.StrategyMode }
    {Qu(OpState)		PT.Control.OperationState }
    {Qu(FuelLevel)		PT.Engine.FuelLevel }
    {Qu(BatterySOC.HV)		PT.BCU.SOC_HV }
    {Qu(Pwr.Exploited)		PT.Control.PwrExploited }
    {Qu(Time.Global)		Time.Global }
    {Qu(OperationState_Tgt)	DM.OperationState_tgt }
    {Qu(OperationState)		Vhcl.OperationState }
    {Qu(Ignition) 		Vhcl.Ignition }
    {Qu(LapTime)	DM.Lap.Time 	}
    {Qu(LapNo)		DM.Lap.No 	}
    {Qu(TriggerPointTime)	DM.TriggerPoint.Time }
    {Qu(TriggerPointDist)	DM.TriggerPoint.Dist }
    {Qu(TriggerPointId)		DM.TriggerPoint.Id   }

}

set Qu(ACC_IsActive) 0
set Qu(EngineSwitch) 0
set Qu(StratMode.Act) ""
set Qu(OpState.Act) ""
set Qu(Time.Global.Act) 0
set Qu(OperationState.Act) 0
set Qu(Ignition.Act) 0

# max 8 Lights
set LEDs {
    0 "PWR" red
    1 "START" #0C0
    2 "ABS" blue 
    3 "ESP" orange
    4 "ACC" #0C0
    5 "Obs" red
    6 "SPD\n-" #0C0
    7 "SPD\n+" #0C0
}

set LEDs_MM {
    0 "PWR" red   
    2 "ABS" blue 
    4 "(o)" red   
    5 "BLS" red
}



## Handle args of command line ################################################

proc PrintUsage {} {
    global HIL Config
puts "
usage: Instruments \[options\]
 -appclass c     restrict connections to appl class c         ($HIL(AcceptClass))
 -apphost h      .. to appl running on host h n               ($HIL(AcceptHost))
 -appuser u      .. to appl owned by user u                   ($HIL(AcceptUser))
 -apppid pid     .. to appl process pid                       ($HIL(AcceptPid))
 -connect        auto-connect after program start
 -maxspeed v     set max speed                                ($Config(MaxSpeed))
 -maxrpm v       set max RPM                                  ($Config(MaxRPM))
 -maxstw v       set max steering wheel angle                 ($Config(MaxStW))
 -ngears v       set number of gears                          ($Config(nGears))
 -units s        configure system of measurement              ($Config(Units))
 -fueleff e      configure fuel efficiency mode               ($Config(FuelEff))
 -eleceff e      configure energy efficiency mode             ($Config(ElecEff))
 -automatic      select automatic gearbox                     ($Config(Automatic))
 -automanual     select automatic gearbox with manual control ($Config(AutoManual))
 -MM             MotorCycleMaker mode
"
	exit 1
}

for {set i 0} {$i<$argc} {incr i} {
    set arg [lindex $argv $i]
    switch -glob $arg {
     -apphost	 { set HIL(AcceptHost)   [lindex $argv [incr i]] }
     -appuser	 { set HIL(AcceptUser)   [lindex $argv [incr i]] }
     -appclass	 { set HIL(AcceptClass)  [lindex $argv [incr i]] }
     -apppid	 { set HIL(AcceptPid)    [lindex $argv [incr i]] }
     -connect	 { after idle Connect}
     -maxspeed	 { set Config(MaxSpeed)  [lindex $argv [incr i]] }
     -maxrpm	 { set Config(MaxRPM)    [lindex $argv [incr i]] }
     -maxst*	 { set Config(MaxStW)    [lindex $argv [incr i]] }
     -ngears	 { set Config(nGears)    [lindex $argv [incr i]] }
     -units	 { set Config(Units)     [lindex $argv [incr i]] }
     -fueleff    { set Config(FuelEff)   [lindex $argv [incr i]] }
     -eleceff    { set Config(ElecEff) [lindex $argv [incr i]] }
     -automanual { set Config(AutoManual) 1 }
     -automatic	 -
     -auto*	 { set Config(Automatic) 1 }
     -lang	 {  }
     -BM	 -
     -MM         { set Config(Mode)     MM 
     		   set Config(MaxSpeed) 240
		   set Config(MaxRPM)   12000
		   set Config(MaxStW)   50
		   set Config(nGears)   0
		   set LEDs $LEDs_MM
     	         }
     -TM         { set Config(Mode)     TM
     		   set Config(MaxSpeed) 125
		   set Config(MaxRPM)   3000
		   set Config(MaxStW)   540
		   set Config(nGears)   0
     	         }
     -*          { PrintUsage }
    }
}

set Config(RPMscale) [expr {$Config(MaxRPM) > 5500 ? 1000 : 100}]

proc MM {a b} {
    return [expr {$::Config(Mode)=="MM" ? $a : $b}]
}


## Load helper modules #######

if {![file exists $Pgm(Dir)/Tacho.tcl]} { append Pgm(Dir) "/../GUI/lib" } 
if {![file exists $Pgm(Dir)/Tacho.tcl]} { set Pgm(Dir) "$env(CM_HOME)/GUI/lib" }

load $LibDir/apoc.so

source $Pgm(Dir)/Tacho.tcl
source $Pgm(Dir)/HVBar.tcl
source $Pgm(Dir)/Gearbox.tcl
source $Pgm(Dir)/MultiPanel.tcl


### Communication port for remote control (send/dde) ###

tk appname Instruments
if {$Pgm(Arch)=="win32"} { dde servername Instruments }


## APO communication ##################################################################

set ApoId ""


proc Subscribe {{sid ""}} {
    global ApoId HIL Qu Subscription Missing

    foreach qu [dataprocess dictionary $ApoId] {
    	set Dict([lindex $qu 0]) $qu
    }
    foreach sub $Subscription {
	set var [lindex $sub 0]
	for {set i 1} {$i < [llength $sub]} {incr i} {
	    set qu [lindex $sub $i]
	    if {[info exist Dict($qu)]} break
	}
        lappend Subs [list $qu $var]
	set $var 0
    }
    set Missing {}
    set Rejected [dataprocess subscribe $ApoId 0 100 HandleDataVec $Subs]
    foreach i $Rejected {
	puts "Not available: #$i [lindex $Subs $i 0]"
	lappend Missing [lindex $Subs $i 1]
    }
}


proc Connect {{Target ""}} {
    global ApoId Pgm HIL

    if {$Target != ""} {
	set HIL(AcceptPid) ""
	if {[scan $Target {%[^: ]:%d} HIL(AcceptHost) HIL(AcceptPid)] != 2} {
	    if {[scan $Target {%u.%u} HIL(AcceptPid) dummy] == 1} { 
	    	set HIL(AcceptHost) localhost
	    }
	}
    }

    apoc query 1000 $HIL(AcceptHost)
    while {![apoc querydone]} {after 100; update}

    set SelIndex      -1
    set Sel(AppStart) 0
    set Sel(AppHost)  ""
    set i 0
    foreach s [apoc servers] {
	lassign $s Info(AppClass) Info(AppName) Info(AppInfo) dummy Info(AppDescr) \
		   Info(AppStart) Info(AppUser) Info(AppHost) Info(AppPid)
	if {($HIL(AcceptUser)==""  || $HIL(AcceptUser)==$Info(AppUser))
	 && ($HIL(AcceptClass)=="" || $HIL(AcceptClass)==$Info(AppClass))
	 && ($HIL(AcceptPid)==""   || $HIL(AcceptPid)==$Info(AppPid))
	 && $Info(AppStart) > $Sel(AppStart)} {
	    set SelIndex $i
	    array set Sel [array get Info]
	}
	incr i
    }
    if {$SelIndex < 0} {return 0}
    Disconnect

    if {[catch {dataprocess new $SelIndex 33} ApoId] != 0} {
	wm title . "$Pgm(Name) - could not connect"
	set ApoId ""
	return
    }
    dataprocess tracemsg $ApoId 5 RecvMsg
    dataprocess dictchanged $ApoId Subscribe
    wm title . "$Pgm(Name) - '$Sel(AppHost)' online"

    Subscribe
    set ::Qu(SpeedLimit) -1
    SendRequest
}


proc Disconnect {} {
    catch {dataprocess delete $::ApoId}
    set ::ApoId ""
    wm title . "$::Pgm(Name) - disconnected"
}


proc SendMsg {Channel Msg} {
    catch {dataprocess sendmsg $::ApoId $Channel $Msg}
}


proc RecvMsg {Process ch Msg} {
    binary scan $Msg S Id
    if {$Id == 0} {
     	set n [binary scan $Msg x2ccSSRRSS AutoGBUsed BrakeKind nFGears nBGears totvEngineMax velVhclMax PTKind GearBoxKind]
	set ptkind [expr {$n==6 ? 0 : $PTKind}]
	set auto   [expr {$n==6 ? $AutoGBUsed : ($GearBoxKind != 1)}]
	set manual [expr {$n==6 ? $AutoGBUsed : ($GearBoxKind == 2)}]
	Gearbox::Configure .gear.no -n $nFGears -auto $auto -automanual $manual
	Tacho::Configure .rpm -max $totvEngineMax
	Tacho::Configure .spd -max $velVhclMax
	HVBar::Configure .mid.clutch -auto $auto
	Configure $ptkind
    }
}


proc SendDVAwr {QuName Value {Where 0} {Mode 0} {nActive 1}} {
    set msg [binary f a8I "DVAwr   " 123]
    # Type(f/i) Where(in/out/DM) Mode(val/offset/factor) nActive QName
    append msg [binary f accIa*xR f $Where $Mode $nActive $QuName $Value]
    SendMsg 2 $msg
}


proc ReleaseDVAwr {{QuName "*"}} {
    set msg  [binary f a8I "DVAwr   " 123]
    append msg [binary f accIa*xI f 0 0 0 $QuName 0]
    SendMsg 2 $msg
}


proc SendRequest {} {
    set msg [binary f S 1000]
    SendMsg 5 $msg
}


###############################################################################


proc HandleDataVec {sid msgid} {
    global Config Qu LED Missing

    Tacho::SetNeedle .rpm 0 [expr {9.549*$Qu(nEngine)/$Config(RPMscale)}]
    Tacho::SetNeedle .spd 0 [CalcDist [expr {$Qu(vCar)*3.6}]]
    Tacho::SetNeedle .spd 1 [CalcDist [expr {$Qu(SpeedLimit)*3.6}]]
    Tacho::SetNeedle .pow 0 $Qu(Pwr.Exploited)


    .spd itemconfig km -text [format "%.3f" [CalcDist [expr {$Qu(Distance)*0.001}]]]
    .spd itemconfig indL -fill [expr {$Qu(IndL) > 1 ? {#1e1} : {grey20} }]
    .spd itemconfig indR -fill [expr {$Qu(IndR) > 1 ? {#1e1} : {grey20} }]
    .rpm itemconfig gear -text [Gearbox::GearName $Qu(PT.GearNo)]
    .pow itemconfig gear -text [Gearbox::GearName $Qu(PT.GearNo)]
    .pow itemconfig rpminfo -text [expr {max(0,(int(9.549*$Qu(nEngine))/10)*10)}]

    # consumption display
    .mul itemconfig fuel_act -text [CalcFuelEff $Qu(FuelConsump_Act)]
    .mul itemconfig fuel_avg -text [CalcFuelEff $Qu(FuelConsump_Avg)]
    .mul itemconfig fuel_lvl -text "[expr {int($Qu(FuelLevel))}] %"
    .mul itemconfig elec_act -text [CalcElecEff $Qu(ElecConsump_Act)]
    .mul itemconfig elec_avg -text [CalcElecEff $Qu(ElecConsump_Avg)]
    .mul itemconfig elec_soc -text "[expr {int($Qu(BatterySOC.HV))}] %"

    if {$Qu(vCar)<1} {
	.mul itemconfig fuel_act -text ""
	.mul itemconfig elec_act -text ""
    }
    if {$Qu(Distance)<10} {
	.mul itemconfig fuel_avg -text ""
	.mul itemconfig elec_avg -text ""
    }

    HVBar::Set .mid.steer  [expr {$Qu(StWhlAngle)*57.296}]
    HVBar::Set .mid.gas    [expr {100*$Qu(Gas)}]
    HVBar::Set .mid.brake  [expr {100*$Qu(Brake)}]
    HVBar::Set .mid.brake2 [expr {100*$Qu(BrakeLever)}]
    HVBar::Set .mid.clutch [expr {100*$Qu(Clutch)}]
    if {$Config(Automatic)} {
	Gearbox::Set .gear.no $Qu(SelectorCtrl)
    } elseif {$Config(nGears)>0} {
	Gearbox::Set .gear.no $Qu(GearNo)
    }

    if {$Qu(Ignition) != $Qu(Ignition.Act)} {
	.lights.l0.l configure -bg $LED(0,$Qu(Ignition))
	set Qu(Ignition.Act) $Qu(Ignition)
    }
    if {$Qu(OperationState) != $Qu(OperationState.Act)} {
	.lights.l1.l configure -bg $LED(1,[expr {$Qu(OperationState) == 4}])
	set Qu(OperationState.Act) $Qu(OperationState)
    }

    .lights itemconfigure TrapTime -text [format " %.3f s"       $Qu(TriggerPointTime)]
    .lights itemconfigure TrapId   -text [format " Sector:   %i" $Qu(TriggerPointId)]
    .lights itemconfigure LapTime  -text [format " %.3f s"       $Qu(LapTime)]
    .lights itemconfigure LapNo    -text [format " Lap Time: %i" $Qu(LapNo)]


    if {$Qu(StratMode) != $Qu(StratMode.Act) || $Qu(OpState.Act) != $Qu(OpState)} {
	MultiPanel::AdjustPTFlow .mul $Qu(StratMode) $Qu(OpState)
	set Qu(StratMode.Act) $Qu(StratMode)
	set Qu(OpState.Act) $Qu(OpState)
    }
    set Qu(UpdatePTFlow) [expr {$Qu(Time.Global.Act)==$Qu(Time.Global) || "Qu(StratMode)" in $Missing ? 0 : 1}]
    set Qu(Time.Global.Act) $Qu(Time.Global)
}


################################################################################


proc CreateDialog {} {
    global simData Config Pref Pgm
    upvar #0 simData v

    if {$Config(FuelEff) eq ""} { set Config(FuelEff) $Pref(FuelEff) }
    if {$Config(ElecEff) eq ""} { set Config(ElecEff) $Pref(ElecEff) }
    if {$Config(Units)   eq ""} { set Config(Units) $Pref(Units) }

    TopLevel . "$Pgm(Name) $Pgm(Version) $Pgm(VersionDate)"
    wm iconname . $Pgm(Name)
    wm protocol  . WM_DELETE_WINDOW ProgramExit
    wm resizable . 0 0

    frame .top -relief raised -bd 1

    ## Tacho
    Tacho::Create .spd -width 220 -height 190 -scaledist 1.17 \
	    -max $Config(MaxSpeed) -majorlen 0.1 -minorlen 0.05 \
	    -major [expr {$Config(MaxSpeed)<130 ? 10 : 20}] -nminor 4 \
            -cap 0 -pad 28 -font "Helvetica -15 bold" -cpad black \
            -bar 0 {0.78  0.82 grey20 grey20} \
            -needle 0 {0.78 1.10 4 red} \
            -needle 1 {0.94 0.99 3 red} \
            -needle 2 {0.88 0.92 5 #0C0} 

    .spd create text 110 65 -text "km/h" -anchor center -fill white \
	    -font {Sans -12} -tags {bg speed}
    .spd create rect 82 169 139 186 -fill grey20 -outline grey40 -tags bg
    .spd create text 136 178 -text "0.000" -anchor e -fill grey90 \
	    -font {Sans -14 bold} -tags km
    .spd create rect 70 95 77 105 -fill grey20 -width 0 -tags indL
    .spd create polygon 60 100 70 90 70 110 -fill grey20 -width 0 -tag indL 
    .spd create rect 143 95 150 105 -fill grey20 -width 0 -tags indR
    .spd create polygon 160 100 150 90 150 110 -fill grey20 -width 0 -tag indR 
    .spd lower bg needles

    ## Multi Panel
    MultiPanel::Create .mul -width 160 -height 190 

    ## Motor RPM
    Tacho::Create .rpm -width 220 -height 190 -scaledist 1.17 \
	    -max   [expr {$Config(MaxRPM)/$Config(RPMscale)}] \
	    -major [expr {$Config(MaxRPM)<4000 ? 5 : 1000/$Config(RPMscale)}] \
	    -nminor 5 -cap 0 -pad 28 -font "Helvetica -15 bold" -cpad black \
            -bar 0 {0.78  0.82 grey20 grey20} \
            -needle 0 {0.78 1.10 4 red}

    .rpm create text 110 85 -text "[l "rpm" "U/min"]\nx $Config(RPMscale)" \
    	    -anchor center -ju c -fill white -font {Sans -12} -tags {bg rpm}
    .rpm create text 110 147 -text "Gear" -anchor c -fill grey90 -font {Sans -9} -tags bg
    .rpm create rect  95 155 125 186 -fill grey20 -outline grey40 -tags bg
    .rpm create text 110 172 -anchor c -fill grey90 -font {Sans -24 bold} -tags {bg gear}
    .rpm lower bg needles

    ## Powermeter
    Tacho::Create .pow -width 220 -height 190 -scaledist 1.17 \
	    -min -20 -max 120 -major 20 -nminor 4 \
            -cap 0 -pad 28 -font "Helvetica -15 bold" -cpad black \
            -bar 0 {0.78  0.82 grey20 grey20} \
            -needle 0 {0.78 1.10 4 red}

    .pow create text 110 65 -text "% PWR" \
    	    -anchor center -ju c -fill white -font {Sans -12} -tags bg
    .pow create rect  22 163 55 175 -fill black -tags bg
    .pow create rect  165 163 198 175 -fill black -tags bg

    .pow create text 110 120 -text "RPM" -anchor c -fill grey90 -font {Sans -11} -tags {bg rpm_txt}
    .pow create rect  75 90 145 110 -fill grey20 -outline grey40 -tags {bg rpm_txt}
    .pow create text 110 100 -anchor c -fill grey90 -font {Sans -20 bold} -tags {bg rpm_txt rpminfo}

    .pow create text 18 130 -anchor c -fill white -font {Sans -14 bold} -angle -69 -text "C"
    .pow create text 22 140 -anchor c -fill white -font {Sans -14 bold} -angle -66 -text "H"
    .pow create text 25 148 -anchor c -fill white -font {Sans -14 bold} -angle -62 -text "A"
    .pow create text 31 158 -anchor c -fill white -font {Sans -14 bold} -angle -59 -text "R"
    .pow create text 36 168 -anchor c -fill white -font {Sans -14 bold} -angle -56 -text "G"
    .pow create text 42 176 -anchor c -fill white -font {Sans -14 bold} -angle -52 -text "E"

    .pow create text 182 170 -anchor c -fill white -font {Sans -14 bold} -angle 53 -text "B"
    .pow create text 187 162 -anchor c -fill white -font {Sans -14 bold} -angle 56 -text "O"
    .pow create text 193 153 -anchor c -fill white -font {Sans -14 bold} -angle 59 -text "O"
    .pow create text 199 144 -anchor c -fill white -font {Sans -14 bold} -angle 62 -text "S"
    .pow create text 202 137 -anchor c -fill white -font {Sans -14 bold} -angle 65 -text "T"

    .pow create text 110 147 -text "Gear" -anchor c -fill grey90 -font {Sans -9} -tags bg
    .pow create rect  95 155 125 186 -fill grey20 -outline grey40 -tags bg
    .pow create text 110 172 -anchor c -fill grey90 -font {Sans -24 bold} -tags {bg gear}
    .pow lower bg needles

    gridrc 1 1 .spd -in .top -sticky w
    gridrc 1 3 .mul -in .top
    gridrc 1 5 .rpm -in .top -sticky e

    grid rowconfigure .top 0 -minsize 5
    grid columnconfigure .top 0 -minsize 3
    grid columnconfigure .top 2 -minsize 3
    grid columnconfigure .top 4 -minsize 3
    grid columnconfigure .top 6 -minsize 3

    ## Control Lights
    global LED LEDs
    set cv .lights
    canvas $cv -width 120 -height 170 -relief raised -bd 1
    foreach {i txt co} $LEDs {
	set LED($i,0)   grey70
	set LED($i,1)   $co
	set LED($i,txt) $txt
    }
    for {set i 0} {$i<2} {incr i} {
	if {![info exists LED($i,0)]} {
	    set LED($i,0)   grey70
	    set LED($i,1)   red
	    set LED($i,txt) ""
	}
	set x [expr {15 + $i%2 * 47}]
	set y [expr {9 + $i/2 * 40}]
	frame $cv.l$i -relief sunken -bd 1
	label $cv.l$i.l -fg white -text $LED($i,txt) -font {Sans -12 bold}
	pack  $cv.l$i.l -fill both -expand 1
	bind  $cv.l$i.l <ButtonRelease> {%W configure -relief flat}
	$cv create window $x $y -anchor nw -width 42 -height 35 -window $cv.l$i
    }

    bind $cv.l0.l <ButtonPress>   { if {![info exists Qu(OperationState)]} return
				    SendDVAwr DM.OperationState_trg [expr {$Qu(Ignition) ? 1 : 3}];
				    %W configure -relief sunken -bg #fcc
				  }
    bind $cv.l1.l <ButtonPress>   { if {![info exists Qu(OperationState)]} return
				    SendDVAwr DM.OperationState_trg [expr {$Qu(OperationState) != 4 ? 4 : 3}];
				    %W configure -relief sunken -bg #cec
				  }
    
	## Laptime
        .lights create text 60  80 -text 0 -anchor center -fill blue  -font {Helvetica -16 bold} -tags TrapId
        .lights create text 60 100 -text 0 -anchor center -fill black -font {Helvetica -14 bold} -tags TrapTime
        .lights create text 60 120 -text 0 -anchor center -fill blue  -font {Helvetica -16 bold} -tags LapNo
        .lights create text 60 140 -text 0 -anchor center -fill red   -font {Helvetica -14 bold} -tags LapTime
	

    frame .mid -relief raised -bd 1

    ## Steering Wheel
    set txt [MM [l "Steering Angle \[deg]" "Lenkwinkel \[Grad]"] \
		[l "Steering Wheel Angle \[deg]" "Lenkradwinkel \[Grad]"] ]
    HVBar::Create .mid.steer -mode needle -width [MM 280 250] -height 64 \
	    -title $txt -color #1e1 \
	    -min -$Config(MaxStW) -max $Config(MaxStW) -inverse \
	    -major [MM 10 90] -nminor [MM 5 3]
    
    ## Clutch
    HVBar::CreateVBar .mid.clutch -width 65 -height 100 -major 20 -color #11f \
	    -title [l "Clutch" "Kupplung"] -unit "%" 
    
    ## Gas Pedal
    HVBar::CreateVBar .mid.gas -width 65 -height 100 -major 20 -color #1e1 \
	    -title "Gas" -unit "%" 
    
    ## Brake Pedal
    HVBar::CreateVBar .mid.brake -width 75 -height 100 -major 20 -color #f11 \
	    -title [MM [l "Brake Pedal" "Bremspedal"] [l "Brake" "Bremse"]] -unit "%"
    
    ## Brake Lever
    HVBar::CreateVBar .mid.brake2 -width 75 -height 100 -major 20 -color #f11 \
	    -title [l "Brake Lever" "Bremshebel"] -unit "%"
    
    ## Gearbox
    frame .gear -relief raised -bd 1
    if {$Config(Mode)!="MM" && $Config(nGears)>0} {
	label .gear.la -font SansBold -text [l "Gearbox" "Schaltung"] 
	Gearbox::Create .gear.no -width 120 -height 142 -n $Config(nGears) \
		    -auto $Config(Automatic) -automanual $Config(AutoManual)
	pack .gear.la -side top -anchor n -padx 3
	pack .gear.no -side top -anchor n -expand 1
    }

    grid .top -columnspan 6 -sticky we

    grid .lights .mid .gear -sticky news
    grid columnconfig . 1 -weight 1
    grid columnconfig . 2 -minsize 100 -weight 1

    grid x .mid.steer - - - -sticky news
    if {$Config(Mode)=="MM"} {
	grid x .mid.clutch .mid.brake .mid.brake2 .mid.gas 
    } {
	grid x .mid.clutch .mid.brake .mid.gas -sticky ews
    }
    grid columnconfig .mid 0 -weight 1 
    grid columnconfig .mid 5 -weight 1

    ## Initialize Instruments
    Tacho::SetNeedle .spd 0 0
    Tacho::SetNeedle .rpm 0 0
    Tacho::SetNeedle .pow 0 0

    ## Context menu
    set m [menu .cmenu -tearoff 0]
    $m add radiobutton -lab "Metric units" -variable Config(Units) \
	-value "metric" -command Configure
    $m add radiobutton -lab "Imperial units" -variable Config(Units) \
	-value "imperial" -command Configure
    $m add radiobutton -lab "US-American units" -variable Config(Units) \
	-value "american" -command Configure
    $m add separator
    $m add cascade -lab "Fuel Efficiency" -menu $m.fuel
    menu $m.fuel -tearoff 0
    $m.fuel add radiobutton -lab "litres / 100 km" -variable Config(FuelEff) \
	-value "consumption" -command Configure
    $m.fuel add radiobutton -lab "km / litre" -variable Config(FuelEff) \
	-value "economy" -command Configure
    $m add cascade -lab "Energy Efficiency" -menu $m.energy
    menu $m.energy -tearoff 0
    $m.energy add radiobutton -lab "kWh / 100 km" -variable Config(ElecEff) \
	-value "consumption" -command Configure
    $m.energy add radiobutton -lab "km / kWh" -variable Config(ElecEff) \
	-value "economy" -command Configure


    $m add separator
    $m add com -lab "Connect"    -u 0 -comm Connect
    $m add com -lab "Disconnect" -u 0 -comm Disconnect
    $m add com -lab "Exit"       -u 2 -comm ProgramExit
    bind all <Button-3> "tk_popup .cmenu %X %Y"
    Configure
}


proc Configure {{PTKind ""}} {
    global Config

    switch $PTKind {
      0 { set Config(PTKind) "Generic" }
      1 { set Config(PTKind) "HEV_Parallel_P1" }
      2 { set Config(PTKind) "HEV_Parallel_P2" }
      3 { set Config(PTKind) "HEV_AxleSplit" }
      4 { set Config(PTKind) "HEV_PowerSplit" }
      5 { set Config(PTKind) "HEV_Serial" }
      6 { set Config(PTKind) "BEV" }
    }

    Tacho::Configure .spd -max $Config(MaxSpeed)

    if {$Config(Units) eq "metric"} {
	set stext "km/h"
	if {$Config(FuelEff) eq "consumption"} {
	    set ltext "l / 100 km"
	} else {
	    set ltext "km / l"
	}
	.cmenu.fuel entryconfigure 0 -lab "litres / 100 km"
	.cmenu.fuel entryconfigure 1 -lab "km / litre"
	if {$Config(ElecEff) eq "consumption"} {
	    set ltext2 "kWh / 100 km"
	} else {
	    set ltext2 "km / kWh"
	}
	.cmenu.energy entryconfigure 0 -lab "kW / 100 km"
	.cmenu.energy entryconfigure 1 -lab "km / kW"
    } else {
	set stext "MPH"
	if {$Config(FuelEff) eq "consumption"} {
	    set ltext "GPM"
	} else {
	    set ltext "MPG"
	}
	.cmenu.fuel entryconfigure 0 -lab "gallons / 100 miles"
	.cmenu.fuel entryconfigure 1 -lab "miles / gallon"
	if {$Config(ElecEff) eq "consumption"} {
	    set ltext2 "kWh / 100 M"
	} else {
	    set ltext2 "M / kWh"
	}
	.cmenu.energy entryconfigure 0 -lab "kWh / 100 miles"
	.cmenu.energy entryconfigure 1 -lab "miles / kWh"
    }
    .mul itemconfigure fuel_txt -text $ltext
    .mul itemconfigure elec_txt -text $ltext2
    .spd itemconfigure speed -text $stext

    if {$Config(PTKind) in { "Generic" }} {
	grid remove .pow
	gridrc 1 5 .rpm -in .top -sticky e
    } else {
	grid remove .rpm
	gridrc 1 5 .pow -in .top -sticky e
    }

    if {$Config(PTKind) in { "HEV_Serial" "BEV" }} {
	.pow itemconfigure rpm_txt -state hidden
    } else {
	.pow itemconfigure rpm_txt -state normal
    }

    MultiPanel::SetPTKind .mul $Config(PTKind)
}


proc CalcFuelEff {val} {
    global Config
    if {$Config(FuelEff) == "economy" && $val == 0} { return "" }
    switch [list $Config(Units) $Config(FuelEff)] {
	{metric consumption}	{ }
	{metric economy}	{ set val [expr {100.0/$val} ] } 
	{imperial consumption}	{ set val [expr {$val/2.8248} ] }
	{imperial economy}	{ set val [expr {282.48/$val} ] }
	{american consumption}	{ set val [expr {$val/2.3521} ] }
	{american economy}	{ set val [expr {235.21/$val} ] }
    }
    return [format "%.2f" [expr {min(99.99, $val)} ]]
}


proc CalcElecEff {val} {
    global Config
    if {$Config(ElecEff) == "economy" && $val == 0} { return "" }
    switch [list $Config(Units) $Config(ElecEff)] {
	{metric consumption}	{ }
	{metric economy}	{ set val [expr {100.0/$val} ] } 
	{imperial consumption}	-
	{american consumption}	{ set val [expr {$val*1.609} ] }
	{imperial economy}	-
	{american economy}	{ set val [expr {62.1504/$val} ] }
    }
    return [format "%.2f" [expr {max(0.0, min(99.99, $val))} ]]
}


proc CalcDist {val} {
    return [expr {$::Config(Units) eq "metric" ? $val : 0.621*$val} ]
}


## Preferences ###################################

array set Pref {			
    Window.Main	""
    FuelEff     "consumption"
    ElecEff     "consumption"
    Units       "metric"
}
catch {source ~/.instruments.rc}


proc PrefSave {} {
    global Pref Config
    scan [wm geo .] "%*dx%*d%s" Pref(Window.Main)
    set Pref(FuelEff) $Config(FuelEff)
    set Pref(Units) $Config(Units)
    catch {
        set fd [open "~/.instruments.rc" w]
        foreach n [array names Pref] { puts $fd "set Pref($n) $Pref($n)" }
	close $fd
    }
}


## MAIN ##########################################

proc ProgramExit {} {
    PrefSave
    exit
}

CreateDialog
bind all <Key-F2> 	{SendMsg 2 "Kl15\0"}
bind all <Control-Key-r> Connect
bind all <Control-Key-q> ProgramExit
bind all <Alt-Key-x>     ProgramExit

apoc pollinterval 100
