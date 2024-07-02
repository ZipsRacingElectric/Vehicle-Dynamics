##
## IndicaterLights.tcl
## CarMaker 7.0 ScriptControl Example - IPG Automotive GmbH (www.ipg-automotive.com)
##
## Example showing how to time the Indicator/Fault Lamps
## to determine if an error state exists
##
## REMARK:
## For this script to run properly, IO quantities need to be provided by the HIL backend.
##
## $Id: IndicatorLights.tcl,v 1.5 2017/04/11 09:12:57 vw Exp $

Log "\n** Indicator Light Test **\n"

# Subscribe to the Lamp quantities
QuantSubscribe { Time Vhcl.Ignition IO.LaABS IO.LaASRFDR IO.LaBrake }
Sleep 3000


# Create the clocks that will be used to time the lamps
ClockCreate ABSClock {$Qu(IO.LaABS) > 0}
ClockCreate ESPClock {$Qu(IO.LaASRFDR) > 0}
#ClockCreate EBVClock {$Qu(IO.LaBrake) > 0}


# Turn the power on.
# If the power is already on, turn it off, wait two seconds,
# and turn it on again.

if {$Qu(Vhcl.Ignition) == 0} {
    PowerOn
} else { 
    PowerOff
    Sleep 2000
    PowerOn
}


# Wait for the lights to turn on before starting
# the clock timers

Log "Waiting for Lights to turn on ..."
WaitForCondition {$Qu(IO.LaABS)    > 0}   10000
WaitForCondition {$Qu(IO.LaASRFDR) > 0}   10000

ClockStart ABSClock
ClockStart ESPClock

# While the lights are on, stay in the loop.
# When one of the lights turns off, stop its clock.
# When all the lights are off, drop out of loop

Log "Timing duration ..."
while {$Qu(IO.LaABS) || $Qu(IO.LaASRFDR) } { 
    if {!$Qu(IO.LaABS)} {
        ClockStop ABSClock
    }
    if {!$Qu(IO.LaASRFDR)} {
        ClockStop ESPClock
    }
    Sleep 5
}


# Get the times of each light

set ABSDuration [ClockGetTime ABSClock]
set ESPDuration [ClockGetTime ESPClock]

# Set the maximum duration a light should stay on
# before an error condition is assumed
set MAXDuration 3


# Check to see if the light was on longer than
# the value of MAXDuration. Print the time 
# and a "Good" or "Bad" message

Log "++++++++++++++++++++++++++++++++++++++++++++++++"

if {$ABSDuration <= $MAXDuration} {
    Log "    Duration: $ABSDuration = ABS GOOD"
} else {
    Log "    Duration: $ABSDuration = ABS BAD"
}

if {$ESPDuration <= $MAXDuration} {
    Log "    Duration: $ESPDuration = ESP GOOD"
} else {
    Log "    Duration: $ESPDuration = ESP BAD"
}


# don't forget to free resources after use
ClockDelete ABSClock
ClockDelete ESPClock


Log "\n** END Indicator Light Test **\n"

PowerOff
