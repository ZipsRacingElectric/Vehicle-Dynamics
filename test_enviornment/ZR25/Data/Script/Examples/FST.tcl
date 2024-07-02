##
## FST.tcl
## CarMaker 7.0 ScriptControl Example - IPG Automotive GmbH (www.ipg-automotive.com)
##
## This example illustrates the use of the Fail Safe Tester (FST) commands. 
##  
## In Verbose mode, you can see the output in the CarMaker log file.
##  
## The Sleep function are used to slow down execution
## and make it more viewable. They are normally not necessary.
##
## $Id: FST.tcl,v 1.4 2017/04/11 09:12:57 vw Exp $

Log "* Start FST Example\n"


# initialize the FST
FST_Init
Log "Initializing FST ..."
Sleep 1000


# add a channel to a group
#   argument 1 = group   (normally G0, G1, G2 )
#   argument 2 = channel (from the ECUParameters file)   
FST_ChAdd G0 Ignition
Log "Adding Channel Ignition to group G0..."
Sleep 1000


# cut a channel
#   argument 1 = channel (channel must be in a group)
FST_ChCut Ignition
Log "Cutting ignition..."
Sleep 1000


# connect a channel
#  argument 1 = channel (channel should be cut )
FST_ChCon Ignition
Log "Reconnecting ignition ..."
Sleep 1000


# add signal connector
#   argument 1 = signal connector
#   argument 2 = group 
FST_SCAdd 1 G0
Log "Connecting signal connector..."
Sleep 1000


# remove a channel from a group
#   argument 1 = Group
#   argument 2 = Channel
FST_ChRem G0 Ignition
Log "Removing Ignition from group G0..."
Sleep 1000


# add another channel to the group
#   argument 1 = group   (should be G0, G1, G2)
#   argument 2 = channel (from the ECUParameters file)   
FST_ChAdd G0 BLS
Log "Adding Channel BLS to group G0..."
Sleep 1000


# remove a group
#   argument 1 = Group
FST_GrpRem G0
Log "Removing group G0..." 
Sleep 1000 

# reset the FST
FST_Reset
Log "Reseting the FST..."


Log "\n*** END FST.tcl Example *** \n"

