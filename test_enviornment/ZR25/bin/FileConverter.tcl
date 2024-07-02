#! /bin/sh
# next line restarts using wish \
exec CM -wish "$0" -- "$@"

# 2012 IPG Karlsruhe, Josef Henning, Christoph Schick
# File Converter for IPG result files: .erg/.dat --> .dat/.erg
# File Converter for tire data files: .tdx --> IPG tire


# search for GUI "bin"-directory (current/below/above)
if {[file exists FC_Resutil.tcl]} {
    set gui_dir "[pwd]"
} elseif {[file exists "bin/FC_Resutil.tcl"]} {
    set gui_dir "[pwd]/bin"
} else {
     puts stderr "Error!\nCannot find FC_Resutil.tcl!\n"
     tk_messageBox -message "Error!\nCannot find FC_Resutil.tcl!\n"
     exit
}


source $gui_dir/FC_Resutil.tcl
source $gui_dir/FC_Tireutil.tcl

namespace eval FC {
    variable tw .fc

### draw tw GUI
proc Popup {} {
    variable tw
    variable Fc

    if {[winfo exists $tw]} {
        catch {wm deiconify $tw}
        raise $tw
        return
    }

    # title the GUI
    set Title [l "File Converter" "Datei Konvertierer"]
    TopLevel $tw $Title

    wm minsize $tw 390 195
    # no dynamic window adjustment possible
    wm resizable $tw 0 0
    # close the toplevel window
    wm protocol $tw WM_DELETE_WINDOW FC::Close

    # label the GUI and add close button
    label $tw.title -text $Title -font Title -fg darkblue
    button $tw.close -text [l " Close " "Schlieﬂen"] -comm FC::Close

    # define the type frame to choose the conversion format
    FrameWithTitle $tw.fwhat [l "Select type" "W‰hle Format"]
    # destroyOptions: check when the file format changes whether the option GUI should be destroyed
    radiobutton $tw.fwhat.res -text [l "Result files   \[ .erg / .dat\ ]" \
                                    "Ergebnisdateien \[ .erg / .dat\ ]"] \
                              -value "res" -var FC::Fc(Type) -command FC::destroyOptions
    radiobutton $tw.fwhat.tire -text [l "Tire files   \[ .tdx\ ]" "Reifendateien \[ .tdx\ ]"] \
                               -value "tire" -var FC::Fc(Type) -command FC::destroyOptions

    gridrc 1 1 $tw.fwhat.res -sticky w -padx {20 0} -pady {9 3}
    gridrc 1 2 $tw.fwhat.tire -sticky e -padx {0 30} -pady {9 3}
    grid columnconfigure $tw.fwhat 2 -weight 1

    # define the Input/Output frame to load/name files for the GUI
    FrameWithTitle $tw.fio [l "Input / Output" "Eingabe / Ausgabe"]
    label $tw.fio.lin -text [l "Input file " "Eingabedatei"]
    label $tw.fio.lout -text [l "Output file" "Ausgabedatei"]
    entry $tw.fio.in -textvariable FC::Fc(InFile)
    entry $tw.fio.out -textvariable FC::Fc(OutFile)
    button $tw.fio.b_file -relief raised -bd 1 -image [IMG folder_open.16.png] \
                          -width 22 -height 14  -command "FC::BrowseFile {open}"
    button $tw.fio.s_file -relief raised -bd 1 -image [IMG folder_open.16.png] \
                          -width 22 -height 14  -command "FC::BrowseFile {save}"

    gridrc 1 1 $tw.fio.lin  -sticky w -padx 3 -pady {9 0}
    gridrc 1 2 $tw.fio.in -sticky we -pady {9 0}
    gridrc 1 3 $tw.fio.b_file -padx 3 -pady {9 0}
    gridrc 2 1 $tw.fio.lout -sticky w -padx 3 -pady {0 3}
    gridrc 2 2 $tw.fio.out -sticky we -pady {0 3}
    gridrc 2 3 $tw.fio.s_file -padx 3
    grid columnconfigure $tw.fio 2 -weight 1

    # add a checkbutton to open/close the option GUI and add a button to start the conversion
    checkbutton $tw.additional -text [l " Additional options" " Weitere Einstellungen"] \
                               -comm FC::CreateAddOpt
    button $tw.convert -text [l "Converting" "Konvertieren"] -comm FC::Convert -bg #aea

    # draw the tw GUI (without options)
    gridrc 1 1 $tw.title      -sticky w  -padx 5 -pady {10 15}
    gridrc 1 3 $tw.close      -sticky e  -padx 5
    gridrc 2 1 $tw.fwhat      -sticky we -padx 5 -columnspan 3
    gridrc 3 1 $tw.fio        -sticky we -padx 5 -pady 10 -columnspan 3
    gridrc 4 1 $tw.additional -sticky ws  -padx 9 -pady 4
    gridrc 4 3 $tw.convert    -sticky es  -padx 5 -pady {0 5}
    grid columnconfigure $tw 2 -weight 1


}
### end of Popup (draw tw GUI)


proc BrowseFile {OpenSave} {
    variable Fc

    if {$OpenSave == "open"} {
        if {$Fc(Type) == "res"} {
            OpenResFile
        } elseif {$Fc(Type) == "tire"} {
            OpenTireFile
        }
    } elseif {$OpenSave == "save"} {
        if {$Fc(Type) == "res"} {
            SaveResFile
        } elseif {$Fc(Type) == "tire"} {
            SaveTireFile
        }
    }
}


### draw option tw GUI
proc CreateAddOpt {} {
    variable tw
    variable Fc
    # Additional Options active?
    # which kind of option GUI sould be drawn (result/tire)?
    if {$Fc(Type) == "res"} {
        # draw, hide or destroy result option GUI?
        if {$Fc(options)=="dest"} {          ;# option GUI is currently destroyed!
            CreateResOptions                 ;# draw a new result option GUI
            SelAllQuants                     ;# and select all listbox entries
            set Fc(options) "draw"           ;# set option GUI is drawn

        } elseif {$Fc(options)=="draw"} {    ;# option GUI is currently drawn!
            grid remove $tw.fr_opt           ;# remove (hide) result option - Time Base
            grid remove $tw.fr_var           ;# remove (hide) result option - OutputQuantities
            set Fc(options) "hidd"           ;# set option GUI is hidden

        } elseif {$Fc(options)=="hidd"} {    ;# option GUI is currently hidden!
            gridrc 5 1 $tw.fr_opt            ;# draw result option GUI - Time Base
            gridrc 6 1 $tw.fr_var            ;# draw result option GUI - OutputQuantities
            set Fc(options) "draw"           ;# set option GUI is drawn
        }

    } elseif {$Fc(Type)=="tire"} {
        # draw or destroy tire option GUI?
        if {$Fc(options)=="dest"} {
            CreateTireOptions                ;# draw a new tire option GUI
            set Fc(options) "draw"

        } elseif {$Fc(options)=="draw"} {
            destroy $tw.fr_tire              ;# destroy the tire option GUI
            destroy $tw.fr_tire_title
            set Fc(options) "dest"
        }
    }
}
### end of CreateAddOpt (draw option tw GUI)



### destroy the option GUI and clear entry fields when the file type changes
proc destroyOptions {} {
    variable tw
    variable Fc

    # result option GUI existent? / file type changed to tire?
    if {[winfo exists $tw.fr_var]&&$Fc(Type) == "tire"} {
        # destroy/reset all result options (Time Base + OutputQuantities)
        set Fc(InFile) ""
        set Fc(OutFile) ""
        $tw.additional deselect
        ResOptionsReset
        destroy $tw.fr_opt
        destroy $tw.fr_var
        destroy $tw.fr_opt_title
        destroy $tw.fr_var_title
        set Fc(options) "dest"
    # tire option GUI existent? / file type changed to result?
    } elseif {[winfo exists $tw.fr_tire]&&$Fc(Type) == "res"} {
        # destroy/reset all tire options
        set Fc(InFile) ""
        set Fc(OutFile) ""
        $tw.additional deselect
        destroy $tw.fr_tire
        destroy $tw.fr_tire_title
        set Fc(options) "dest"
    }
    # result file loaded? / file type changed to tire?
    if {([file extension $Fc(InFile)] == ".erg"||[file extension $Fc(InFile)] == ".dat") && $Fc(Type) == "tire"} {
        set Fc(InFile) ""
        set Fc(OutFile) ""
        ResOptionsReset
    # tire file loaded? / fite type changed to result?
    } elseif {[file extension $Fc(InFile)] == ".tdx" && $Fc(Type) == "res"} {
        set Fc(InFile) ""
        set Fc(OutFile) ""
    }
}


### destroy the tw GUI
proc Close {} {
    variable tw
    variable Fc

    if {[winfo exists $tw]} {
        set Fc(InFile) ""
        set Fc(OutFile) ""
        ResOptionsReset
        TireOptionsReset
        $tw.additional deselect
        set Fc(options) "dest"
        SaveWindowPref $tw
        destroy $tw
    }
}


### start result or tire conversion?
proc Convert {} {
    variable Fc

    if {![file exists $Fc(InFile)]} {
        puts stderr [l "Error!\nPlease select input file..." \
                    "Error!\nBitte w‰hle eine Eingabedatei aus..."]
        tk_messageBox -message [l "Error!\nPlease select input file..." \
                    "Error!\nBitte w‰hle eine Eingabedatei aus..."]
        return
    }
    if {$Fc(OutFile)==""} {
        puts stderr [l "Error!\nPlease define output file..." \
                    "Error!\nBitte w‰hle eine Ausgabedatei aus..."]
        tk_messageBox -message [l "Error!\nPlease define output file..." \
                    "Error!\nBitte w‰hle eine Ausgabedatei aus..."]
        return
    }

    if {$Fc(Type) == "res"} {
        convertResult
    } elseif {$Fc(Type) == "tire"} {
        convertTire
    }
}


### main() #################################################################################

### general variables
# path to the Input File
set Fc(InFile) ""
# path to the Output File
set Fc(OutFile) ""
# type of the conversion: "res"=Result, "tire"=Tire
set Fc(Type) "res"
# state flag for the option GUI: "dest"=destroyed, "draw"=drawn, "hidd"=hidden
set Fc(options) "dest"
# path to the CarMaker installation
set Fc(instDir) "$Pgm(Dir)/../"
regsub -all { } $Fc(instDir) {\\ } Fc(instDir)
# path to the Resutil installation
set Fc(pathRes) "$Fc(instDir)/bin/resutil"
# path to the Tireutil installation
set Fc(pathTire) "$Fc(instDir)/bin/tireutil"


### variables for the result options
# for the storage of the time vector (from input file)
set Res(timeBase) ""
# current slider Start position (relativ x-position)
set Res(SliderStart) 10
# current slider End position (relativ x-position)
set Res(SliderEnd) 10
# state flag: which slider was moved last?
set Res(movedLast) "SliderStart"
# state flag: to set initial values once
set Res(initialFlag) 1
# limitation of the time base (minimum)
set Res(limitStart) 0.0
# limitation of the time base (maximum)
set Res(limitEnd) 1000.00
# current value for the start time (calculated from SliderStart)
set Res(valueStart) "$Res(limitStart)"
 # current value for the end time (calculated from SliderEnd)
set Res(valueEnd) "$Res(limitEnd)"
# tmp variable to handle the bind commands
set Res(tmpSlider) "$Res(limitEnd)"
# current slider 3 position
set Res(SliderOffset) 0
# current value for the offset time (equals SliderOffset)
set Res(valueOffset) 0.0
# to set the step size for the time vector
set Res(stepSize) 0.02
# step size from file
set Res(stepSizeFile) 0.02
# for the storage of all channels (shown in listbox)
set Res(allQuants) [l "{Select the channels here!}" "{Selektiere hier die Kan‰le!}"]
# for the storage of the selected channels (highlighted blue)
set Res(selQuants) ""
# to show the amount of the slected quantities (is textvariable!)
set Res(labelQuants) [l "(0 from 0 selected)" "(0 von 0 ausgew‰hlt)"]


### variables for the tire options
# to save the parameters from TDX-file
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
# to adapt the parameters
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

if {[pwd] == "/fs/u3/IPG_FCM/work.chs/FileConverter/Entwicklung_FileConverter/src"} {
    Popup
}


}
### end of namespace FC
