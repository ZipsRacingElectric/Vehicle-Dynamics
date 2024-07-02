# 2012 IPG Karlsruhe, Josef Henning, Christoph Schick

# convertion with tireutil (.tdx --> IPG tire)

namespace eval FC {
    variable tw .fc

### create more tire options
proc CreateTireOptions {} {
    variable tw
    variable Fc
    variable Tr

    FrameWithTitle $tw.fr_tire [l "Tire properties" "Reifen-Eigenschaften"]

    frame $tw.fr_tire.fl
    frame $tw.fr_tire.fr
    frame $tw.fr_tire.fu

    label $tw.fr_tire.fl.l_FricCoef -text [l "Friction coefficient \[-]" "Reibungskoeffizient \[-]"]
    label $tw.fr_tire.fl.l_KinTireR -text [l "Kinematic tire radius \[m]" "Kin. Reifenradius \[m]"]
    label $tw.fr_tire.fl.l_NomTireR -text [l "Nominal tire radius \[m]" "Nom. Reifenradius \[m]"]
    label $tw.fr_tire.fl.l_NomRimR -text [l "Nominal rim radius \[m]" "Nom. Felgenradius \[m]"]
    label $tw.fr_tire.fl.l_NomTireW -text [l "Nominal width \[m]" "Nominelle Breite \[m]"]
    label $tw.fr_tire.fl.l_Stiff -text [l "Vertical stiffness \[N/m]" "Vertikale Steifigk. \[N/m]"]
    label $tw.fr_tire.fl.l_Damp -text [l "Vertical damping \[Ns/m]" "Vert. Dämpfung \[Ns/m]"]

    entry $tw.fr_tire.fl.e_FricCoef -textvariable FC::Tr(FricCoef) -width 7 -justify right
    entry $tw.fr_tire.fl.e_KinTireR -textvariable FC::Tr(KinTireR) -width 7 -justify right
    entry $tw.fr_tire.fl.e_NomTireR -textvariable FC::Tr(NomTireR) -width 7 -justify right
    entry $tw.fr_tire.fl.e_NomRimR -textvariable FC::Tr(NomRimR) -width 7 -justify right
    entry $tw.fr_tire.fl.e_NomTireW -textvariable FC::Tr(NomTireW) -width 7 -justify right
    entry $tw.fr_tire.fl.e_Stiff -textvariable FC::Tr(Stiff) -width 7 -justify right
    entry $tw.fr_tire.fl.e_Damp -textvariable FC::Tr(Damp) -width 7 -justify right

    gridrc 1 1 $tw.fr_tire.fl.l_FricCoef -sticky w -padx {3 0} -pady {9 0}
    gridrc 1 2 $tw.fr_tire.fl.e_FricCoef -sticky e -padx {0 0} -pady {9 0}
    gridrc 2 1 $tw.fr_tire.fl.l_KinTireR -sticky w -padx {3 0} -pady {3 0}
    gridrc 2 2 $tw.fr_tire.fl.e_KinTireR -sticky e -padx {0 0} -pady {3 0}
    gridrc 3 1 $tw.fr_tire.fl.l_NomTireR -sticky w -padx {3 0} -pady {3 0}
    gridrc 3 2 $tw.fr_tire.fl.e_NomTireR -sticky e -padx {0 0} -pady {3 0}
    gridrc 4 1 $tw.fr_tire.fl.l_NomRimR -sticky w -padx {3 0} -pady {3 0}
    gridrc 4 2 $tw.fr_tire.fl.e_NomRimR -sticky e -padx {0 0} -pady {3 0}
    gridrc 5 1 $tw.fr_tire.fl.l_NomTireW -sticky w -padx {3 0} -pady {3 0}
    gridrc 5 2 $tw.fr_tire.fl.e_NomTireW -sticky e -padx {0 0} -pady {3 0}
    gridrc 6 1 $tw.fr_tire.fl.l_Stiff -sticky w -padx {3 0} -pady {3 0}
    gridrc 6 2 $tw.fr_tire.fl.e_Stiff -sticky e -padx {0 0} -pady {3 0}
    gridrc 7 1 $tw.fr_tire.fl.l_Damp -sticky w -padx {3 0} -pady {3 0}
    gridrc 7 2 $tw.fr_tire.fl.e_Damp -sticky e -padx {0 0} -pady {3 0}

    label $tw.fr_tire.fr.l_AspRat -text [l "Aspect ratio \[-]" "Querschnittsverhäl. \[-]"]
    label $tw.fr_tire.fr.l_LongRel -text [l "Longit. relaxation \[m]" "Relaxation längs \[m]"]
    label $tw.fr_tire.fr.l_LatRel -text [l "Lateral relaxation \[m]" "Relaxation quer \[m]"]
    label $tw.fr_tire.fr.l_AligRel -text [l "Aligning relaxat. \[m]" "Relaxat. rückstell \[m]"]
    label $tw.fr_tire.fr.l_RollResK -text [l "Roll resistance kind" "Rollwiderstand Typ"]
    label $tw.fr_tire.fr.l_RollRes -text [l "Roll resistance \[-]" "Rollwiderstand \[-]"]
    label $tw.fr_tire.fr.l_FLoadMax -text [l "Load force max. \[N]" "Maximale Radlast \[N]"]

    entry $tw.fr_tire.fr.e_AspRat -textvariable FC::Tr(AspRat) -width 7 -justify right
    entry $tw.fr_tire.fr.e_LongRel -textvariable FC::Tr(LongRel) -width 7 -justify right
    entry $tw.fr_tire.fr.e_LatRel -textvariable FC::Tr(LatRel) -width 7 -justify right
    entry $tw.fr_tire.fr.e_AligRel -textvariable FC::Tr(AligRel) -width 7 -justify right
    entry $tw.fr_tire.fr.e_RollResK -textvariable FC::Tr(RollResK) -width 7
    entry $tw.fr_tire.fr.e_RollRes -textvariable FC::Tr(RollRes) -width 7 -justify right
    entry $tw.fr_tire.fr.e_FLoadMax -textvariable FC::Tr(FLoadMax) -width 7 -justify right

    gridrc 1 1 $tw.fr_tire.fr.l_AspRat -sticky w -padx {3 0} -pady {9 0}
    gridrc 1 2 $tw.fr_tire.fr.e_AspRat -sticky e -padx {0 3} -pady {9 0}
    gridrc 2 1 $tw.fr_tire.fr.l_LongRel -sticky w -padx {3 0} -pady {3 0}
    gridrc 2 2 $tw.fr_tire.fr.e_LongRel -sticky e -padx {0 3} -pady {3 0}
    gridrc 3 1 $tw.fr_tire.fr.l_LatRel -sticky w -padx {3 0} -pady {3 0}
    gridrc 3 2 $tw.fr_tire.fr.e_LatRel -sticky e -padx {0 3} -pady {3 0}
    gridrc 4 1 $tw.fr_tire.fr.l_AligRel -sticky w -padx {3 0} -pady {3 0}
    gridrc 4 2 $tw.fr_tire.fr.e_AligRel -sticky e -padx {0 3} -pady {3 0}
    gridrc 5 1 $tw.fr_tire.fr.l_RollResK -sticky w -padx {3 0} -pady {3 0}
    gridrc 5 2 $tw.fr_tire.fr.e_RollResK -sticky e -padx {0 3} -pady {3 0}
    gridrc 6 1 $tw.fr_tire.fr.l_RollRes -sticky w -padx {3 0} -pady {3 0}
    gridrc 6 2 $tw.fr_tire.fr.e_RollRes -sticky e -padx {0 3} -pady {3 0}
    gridrc 7 1 $tw.fr_tire.fr.l_FLoadMax -sticky w -padx {3 0} -pady {3 0}
    gridrc 7 2 $tw.fr_tire.fr.e_FLoadMax -sticky e -padx {0 3} -pady {3 0}

    # define right click menu for the more properties field
    label $tw.fr_tire.fu.l_more -text [l "Insert more properties" "Weitere Parameter einfügen"]
    set e [entry $tw.fr_tire.fu.e_more -textvariable FC::Tr(more) -justify left -fg grey]
    set m [menu $tw.fr_tire.fu.m_more -tearoff 0 -takefocus 0]

    $m add command -lab [l "Number of load force values \[-]: -nFLoad <value>"\
                        "Anzahl der Radlastwerte \[-]: -nFLoad <value>"] \
                   -comm "FC::InsertProperties $e -nFLoad"

    $m add command -lab [l "Number of slip values \[-]: -nSlip <value>"\
                        "Anzahl der Schlupfwerte \[-]: -nSlip <value>"] \
                   -comm "FC::InsertProperties $e -nSlip"

    $m add command -lab [l "Number of alpha values \[-]: -nAlpha <value>"\
                        "Anzahl der Alphawerte \[-]: -nAlpha <value>"] \
                   -comm "FC::InsertProperties $e -nAlpha"

    $m add command -lab [l "Number of camber values \[-]: -nCamber <value>"\
                        "Anzahl der Sturzwinkel \[-]: -nCamber <value>"] \
                   -comm "FC::InsertProperties $e -nCamber"

    $m add command -lab [l "Minimum value force load \[N]: -FLoadMin <value>"\
                        "Minimalwert der Radlast \[N]: -FLoadMin <value>"] \
                   -comm "FC::InsertProperties $e -FLoadMin"

    $m add command -lab [l "Minimum value slip \[-]: -SlipMin <value>"\
                        "Minimalwert des Schlupfes \[-]: -SlipMin <value>"] \
                   -comm "FC::InsertProperties $e -SlipMin"

    $m add command -lab [l "Minimum value alpha \[deg]: -AlphaMin <value>"\
                        "Minimalwert von Alpha \[deg]: -AlphaMin <value>"] \
                   -comm "FC::InsertProperties $e -AlphaMin"

    $m add command -lab [l "Minimum value camber \[deg]: -CamberMin <value>"\
                        "Minimalwert des Sturzes \[deg]: -CamberMin <value>"] \
                   -comm "FC::InsertProperties $e -CamberMin"

    $m add command -lab [l "Maximum value force load \[N]: -FLoadMax <value>"\
                        "Maximalwert der Radlast \[N]: -FLoadMax <value>"] \
                   -comm "FC::InsertProperties $e -FLoadMax"

    $m add command -lab [l "Maximum value slip \[-]: -SlipMax <value>"\
                        "Maximalwert des Schlupfes \[-]: -SlipMax <value>"] \
                   -comm "FC::InsertProperties $e -SlipMax"

    $m add command -lab [l "Maximum value alpha \[deg]: -AlphaMax <value>"\
                        "Maximalwert von Alpha \[deg]: -AlphaMax <value>"] \
                   -comm "FC::InsertProperties $e -AlphaMax"

    $m add command -lab [l "Maximum value camber \[deg]: -CamberMax <value>"\
                        "Maximalwert des Sturzes \[deg]: -CamberMax <value>"] \
                   -comm "FC::InsertProperties $e -CamberMax"

    $m add command -lab [l "Vertical offset correction of Fx, Fy, Mz: -OffsetCorr <value>"\
                        "Vertikale Offsetkorrektur von Fx, Fy, Mz: -OffsetCorr <value>"] \
                   -comm "FC::InsertProperties $e -OffsetCorr"

    $m add command -lab [l "Horizontal alpha correction of Fy, Mz: -AlphaCorr <value>"\
                        "Horizontale Alphakorrektur von Fy, Mz: -AlphaCorr <value>"] \
                   -comm "FC::InsertProperties $e -AlphaCorr"

    $m add command -lab [l "Estimated correction of Fy, Mz to origin: -FyOrigCorr <value>"\
                        "Korrektur von Fy, Mz zum geschätzten Ursprung: -FyOrigCorr <value>"] \
                   -comm "FC::InsertProperties $e -FyOrigCorr"

    $m add command -lab [l "Horizontal slip correction of Fx: -SlipCorr <value>"\
                        "Horizontale Schlupfkorrektur von Fx: -SlipCorr <value>"] \
                   -comm "FC::InsertProperties $e -SlipCorr"

    $m add command -lab [l "Estimated correction of Fx to origin: -FxOrigCorr <value>"\
                        "Korrektur von Fx zum geschätzten Ursprung: -FxOrigCorr <value>"] \
                   -comm "FC::InsertProperties $e -FxOrigCorr"

    $m add separator

    $m add command -lab [l "Clear more properties field..." "Lösche weitere Parameter..."] \
                   -comm "FC::InsertProperties $e -clear"
 
   $m add command -lab [l "Show help from Tireutil..." "Zeige die Hilfe zu Tireutil..."] \
                  -comm "FC::ShowHelp"

    button $tw.fr_tire.fu.defPar -text [l "  Load TDX-file parameters  " "  Lade TDX-Dateiparameter  "] \
                                 -comm FC::LoadTDXparam

    gridrc 1 1 $tw.fr_tire.fu.l_more -sticky ws -padx {3 0} -pady {9 0}
    gridrc 2 1 $tw.fr_tire.fu.e_more -sticky ew -padx 3 -pady {0 3}
    gridrc 1 1 $tw.fr_tire.fu.defPar -sticky e -padx {0 3} -pady {3 3}

    gridrc 1 1 $tw.fr_tire.fl -sticky w
    gridrc 1 2 $tw.fr_tire.fr -sticky e
    gridrc 2 1 $tw.fr_tire.fu -sticky we -columnspan 2
    grid columnconfigure $tw.fr_tire.fu 1 -weight 1

    gridrc 5 1 $tw.fr_tire  -sticky we -padx 4 -pady {10 5} -columnspan 3
    grid columnconfigure $tw.fr_tire 2 -weight 1

    # binds to clear and re-configure the more properties field
    bind $e <Button-1> {
        FC::configPropertiesField $FC::tw.fr_tire.fu.e_more
    }

    bind $tw.fr_tire.fu.e_more <Button-3> {
        FC::configPropertiesField $FC::tw.fr_tire.fu.e_more
        tk_popup $FC::tw.fr_tire.fu.m_more %X %Y
    }
}
### end of CreateTireOptions


### open and handle the file browser to select the input tire file
proc OpenTireFile {} {
    variable Fc
    variable Tr

    set types {{{TDX files} .tdx}}
    set tmp_input [tk_getOpenFile -filetypes $types -initialdir "[pwd]/Data/Tire"\
        -title [l "Select tire file..." "Wähle Reifendatei aus..."]]
    # no file selected --> do not change the input file name
    if {$tmp_input!=""} {
        set Fc(InFile) $tmp_input
    } else {
        return
    }
    # set the ouput file name (default = input file name)
    if {$Fc(OutFile) == ""} {
        set Fc(OutFile) [file rootname $Fc(InFile)]
    }
    # scroll to end positions
    .fc.fio.in xview end
    .fc.fio.out xview end
}


### open and handle the file browser to select the output tire file
proc SaveTireFile {} {
    variable Fc

    set types {{{IPG tire files} *}}
    set tmp_output [tk_getSaveFile -filetypes $types -initialfile $Fc(OutFile) -initialdir "[pwd]/Data/Tire"\
        -title [l "Select storage path..." "Wähle Speicherpfad aus..."]]
    # no file selected --> do not change the output file name
    if {$tmp_output!=""} {
        set Fc(OutFile) $tmp_output
    } else {
        return
    }
    .fc.fio.out xview end
}


### load the parameters from the TDX-File by clicking on the button
proc LoadTDXparam {} {
    variable Fc
    variable Tr

    if {[file exists $Fc(InFile)]} {

        # scan the TDX-file with tireutil
        # write the parameters from the shell output into "tmpTireOutput.txt"
        set scanID [open tmpTireOutput.txt w+]

        set cmd $Fc(pathTire)
        lappend cmd "-scan"
        lappend cmd "-if"
        lappend cmd "$Fc(InFile)"
        lappend cmd "-c"
        lappend cmd "-of"
        lappend cmd "tmpTireOutput"

        if {[catch [puts $scanID [exec {*}$cmd]] err]} {
            tk_messageBox -message $err
        }

        close $scanID
        set scanID [open tmpTireOutput.txt]

        # scan the tmp-file and save the needed quantities into Tr(...)
        while {1} {
            set line [gets $scanID]

            # search for occuring errors --> uses the shell output
            if {[regexp {Error} $line]||[regexp {ERROR} $line]} {
                close $scanID
                puts stderr "Error!\nCan not read tire file!\nSee \"tmpTireOutput.txt\"..."
                tk_messageBox -message "Error!\nCan not read tire file!\nSee \"tmpTireOutput.txt\"..."
                return
            }

            set Tr(FricCoef) 1.0
            if [regexp {KinRollRadius} $line] {
                set Tr(KinTireR_File) [lindex $line 2]
                set Tr(KinTireR) $Tr(KinTireR_File)
            }
            if [regexp {NomRadius} $line] {
                set Tr(NomTireR_File) [lindex $line 2]
                set Tr(NomTireR) $Tr(NomTireR_File)
            }
            if [regexp {RimRadius} $line] {
                set Tr(NomRimR_File) [lindex $line 2]
                set Tr(NomRimR) $Tr(NomRimR_File)
            }
            if [regexp {NomWidth} $line] {
                set Tr(NomTireW_File) [lindex $line 2]
                set Tr(NomTireW) $Tr(NomTireW_File)
            }
            if [regexp {Radial.Stiffness} $line] {
                set Tr(Stiff_File) [lindex $line 2]
                set Tr(Stiff) $Tr(Stiff_File)
            }
            if [regexp {Radial.Damping} $line] {
                set Tr(Damp_File) [lindex $line 2]
                set Tr(Damp) $Tr(Damp_File)
            }
            if [regexp {AspectRatio} $line] {
                set Tr(AspRat_File) [lindex $line 2]
                set Tr(AspRat) $Tr(AspRat_File)
            }
            if [regexp {LongFrc.Length} $line] {
                set Tr(LongRel_File) [lindex $line 2]
                set Tr(LongRel) $Tr(LongRel_File)
            }
            if [regexp {SideFrc.Length} $line] {
                set Tr(LatRel_File) [lindex $line 2]
                set Tr(LatRel) $Tr(LatRel_File)
            }
            if [regexp {AlignTrq.Length} $line] {
                set Tr(AligRel_File) [lindex $line 2]
                set Tr(AligRel) $Tr(AligRel_File)
            }
            if [regexp {RollResist.Kind} $line] {
                set Tr(RollResK_File) [lindex $line 2]
                set Tr(RollResK) $Tr(RollResK_File)
            }
            if [regexp {RollResist.Factor} $line] {
                set Tr(RollRes_File) [lindex $line 2]
                set Tr(RollRes) $Tr(RollRes_File)
            }
            if [regexp {FLoadMax} $line] {
                set Tr(FLoadMax_File) [lindex $line 2]
                set Tr(FLoadMax) $Tr(FLoadMax_File)
            }
            if {[eof $scanID]} { 
                close $scanID
                break
            }
        }
        # delete all created tmp-files from scanning
        set tmpDelete [glob tmpTireOutput*]
            for {set i 0} {$i < [llength $tmpDelete]} {incr i 1} {
            file delete [lindex $tmpDelete $i]
        }

    } else {
        puts stderr [l "Error!\nPlease select input file..." \
                    "Error!\nBitte wähle eine Eingabedatei aus..."]
        tk_messageBox -message [l "Error!\nPlease select input file..." \
                               "Error!\nBitte wähle eine Eingabedatei aus..."]
    }
}
### end of LoadTDXparam


### clear the more properties field and set the font from grey to black
proc configPropertiesField {w} {

    # set font of entry to black
    $w configure -fg black
    # initial entry existent?
    if {$FC::Tr(more) == "Rechtsklick" || $FC::Tr(more) == "Right click"} {
        set FC::Tr(more) ""
    }
}


### insert the selected option form the right click menu
proc InsertProperties {w prop} {
    variable Tr

    # clear the more properties entry field?
    if {$prop == "-clear"} {
        set Tr(more) ""
        return
    }
    # filter out the needed characters and append the selected option
    lappend Tr(more) $prop
    lappend Tr(more) "<value>"
    set value [expr [string length $Tr(more)]-7]
    # activate the last appended option
    $w selection range $value end
    $w icursor $value
    focus $w
    $w xview end
}


### under construction --> compromise (chs)
### reset all tire variables 
proc TireOptionsReset {} {
    variable Tr

    # reset all variables
    set Tr(KinTireR_File) ""
    set Tr(NomTireR_File) ""
    set Tr(NomRimR_File) ""
    set Tr(NomTireW_File) ""
    set Tr(Stiff_File) ""
    set Tr(Damp_File) ""
    set Tr(AspRat_File) ""
    set Tr(LongRel_File) ""
    set Tr(LatRel_File) ""
    set Tr(AligRel_File) ""
    set Tr(RollResK_File) ""
    set Tr(RollRes_File) ""
    set Tr(FLoadMax_File) ""
    set Tr(FricCoef) ""
    set Tr(KinTireR) ""
    set Tr(NomTireR) ""
    set Tr(NomRimR) ""
    set Tr(NomTireW) ""
    set Tr(Stiff) ""
    set Tr(Damp) ""
    set Tr(AspRat) ""
    set Tr(LongRel) ""
    set Tr(LatRel) ""
    set Tr(AligRel) ""
    set Tr(RollResK) ""
    set Tr(RollRes) ""
    set Tr(FLoadMax) ""
    set Tr(more) [l "Right click" "Rechtsklick"]
}


### start the conversion with tireutil
proc convertTire {} {
    variable Fc
    variable Tr

    set cmd "$Fc(pathTire)"
    lappend cmd "-convert"
    lappend cmd "-if"
    lappend cmd "$Fc(InFile)"

    lappend cmd "-of"
    lappend cmd "$Fc(OutFile)"

    # write the bin-file into the directory "Mapping"
    file mkdir "[file dirname $Fc(OutFile)]/Mapping"
    lappend cmd "-ofbin"
    lappend cmd "[file dirname $Fc(OutFile)]/Mapping/[file tail $Fc(OutFile)].bin"

    # insert only modified parameters
    if {$Tr(FricCoef) != 1.0 && $Tr(FricCoef) != ""} {
        lappend cmd "-muRoad"
        lappend cmd "$Tr(FricCoef)"
    }
    if {$Tr(KinTireR) != $Tr(KinTireR_File)} {
        lappend cmd "-KinRollRadius"
        lappend cmd "$Tr(KinTireR)"
    }
    if {$Tr(NomTireR) != $Tr(NomTireR_File)} {
        lappend cmd "-NomRadius"
        lappend cmd "$Tr(NomTireR)"
    }
    if {$Tr(NomRimR) != $Tr(NomRimR_File)} {
        lappend cmd "-RimRadius"
        lappend cmd "$Tr(NomRimR)"
    }
    if {$Tr(NomTireW) != $Tr(NomTireW_File)} {
        lappend cmd "-NomWidth"
        lappend cmd "$Tr(NomTireW)"
    }
    if {$Tr(Stiff) != $Tr(Stiff_File)} {
        lappend cmd "-cRadial"
        lappend cmd "$Tr(Stiff)"
    }
    if {$Tr(Damp) != $Tr(Damp_File)} {
        lappend cmd "-kRadial"
        lappend cmd "$Tr(Damp)"
    }
    if {$Tr(AspRat) != $Tr(AspRat_File)} {
        lappend cmd "-AspectRatio"
        lappend cmd "$Tr(AspRat)"
    }
    if {$Tr(LongRel) != $Tr(LongRel_File)} {
        lappend cmd "-Length.LongFrc"
        lappend cmd "$Tr(LongRel)"
    }
    if {$Tr(LatRel) != $Tr(LatRel_File)} {
        lappend cmd "-Length.SideFrc"
        lappend cmd "$Tr(LatRel)"
    }
    if {$Tr(AligRel) != $Tr(AligRel_File)} {
        lappend cmd "-Length.AlignTrq"
        lappend cmd "$Tr(AligRel)"
    }
    # options: TrqLoad, TrqVelocity, FLoad, FVelocity
    if {$Tr(RollResK) != $Tr(RollResK_File) || $Tr(RollRes) != $Tr(RollRes_File)} {
        lappend cmd "-RollResist"
        lappend cmd "$Tr(RollResK)"
        lappend cmd "$Tr(RollRes)"
    }
    if {$Tr(FLoadMax) != $Tr(FLoadMax_File)} {
        lappend cmd "-FLoadMax"
        lappend cmd "$Tr(FLoadMax)"
    }
    if {[winfo exists $FC::tw.fr_tire.fu.e_more] && $Tr(more) != ""} {
        configPropertiesField $FC::tw.fr_tire.fu.e_more
        set cmd [concat $cmd $Tr(more)]
    }

    set scanID [open tmpTireOutput.txt w+]

    # execute tireutil and write the shell output into "tmpTireOutput.txt"
    if {[catch [puts $scanID [exec {*}$cmd]] err]} {
        tk_messageBox -message $err

    } else {
        close $scanID
        set scanID [open tmpTireOutput.txt]

        # search for occuring errors --> uses the shell output
        while {1} {
            set line [gets $scanID]
            if {[regexp {Error} $line]||[regexp {ERROR} $line]} {
                close $scanID
                puts stderr "Error!\nTire file not converted!\nSee \"tmpTireOutput.txt\"..."
                tk_messageBox -message "Error!\nTire file not converted!\nSee \"tmpTireOutput.txt\"..."
                return
            }
            if {[eof $scanID]} { 
                break
            }
        }
        tk_messageBox -message [l "File created successfully!" "Datei erfolgreich konvertiert!"]
    }
    close $scanID
    file delete tmpTireOutput.txt
}
### end of convertTire


### show the help menu by right clicking on the more properties entry field
proc ShowHelp {} {
    variable tw

    # menu already existing?
    if {[winfo exists $tw.help]} {
        catch {wm deiconify $tw.help}
        raise $tw.help
        return
    }

# the help text from tireutil
set helpText "tire data utility vCARMAKER_VER
usage: tireutil \[mode\] \[options\]

 mode:
  -convert            convert tire data
  -scan               2-dimensional tire-scan
  -scangrid           3-dimensional tire-scan

 options:
  -v                  verbose
  -c                  check or list parameters
  -if                 input file
  -of                 output file
  -ofbin              output file 'BinFName' of RTTire

 tire properties:
  -KinRollRadius <x>  kinematic tire roll radius
  -NomRadius <x>      nominal tire radius
  -RimRadius <x>      nominal rim radius
  -NomWidth           nominal tire width
  -AspectRatio <x>    aspect ratio
  -muRoad             road friction coefficient
  -cRadial            vertical stiffness coefficient
  -kRadial            vertical damping coefficient
  -Length.LongFrc     longitudinal relaxation length
  -Length.SideFrc     lateral relaxation length
  -Length.AlignTrq    aligning relaxation length
  -RollResist <s> <val>   roll resistance mode

 conversion and scan parameters:
  -n<q>   <i>         number of values
  -<q>Min <x>         minimum value
  -<q>Max <x>         maximum value
                      <q> := FLoad \[N\], Slip \[-\], Alpha \[deg\], Camber\[deg\]
  -OffsetCorr         vertical offset correction of Fx,Fy,Mz
  -AlphaCorr          horizontal alpha correction of Fy,Mz
  -FyOrigCorr         corrects Fy to the estimated origin and corrects
                      Mz with the equal alpha correction
  -SlipCorr           horizontal slip correction of Fx
  -FxOrigCorr         corrects Fx to the estimated origin
  (default: no force/torque correction)

 scan parameters:
  -FLoad <x>          scan for this FLoad only"

    set Title [l "Help Tireutil" "Tireutil-Hilfe"]
    TopLevel $tw.help $Title
    wm resizable $tw.help 0 0

    text $tw.help.txt -height 45 -width 76
    $tw.help.txt insert end $helpText
    button $tw.help.ok -text " Ok " -default active -command {destroy $FC::tw.help}

    gridrc 1 1 $tw.help.txt
    gridrc 2 1 $tw.help.ok
}
### end of ShowHelp

}
### end of namespace FC


