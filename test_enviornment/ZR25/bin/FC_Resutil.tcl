# 2012 IPG Karlsruhe, Josef Henning, Christoph Schick

# convertion with resutil (convertion of .erg <--> .dat)

namespace eval FC {
    variable tw .fc

### create more result options
proc CreateResOptions {} {
    variable tw
    variable Fc
    variable Res

    FrameWithTitle $tw.fr_opt [l "Time base" "Zeitbasis"]

    # frame 1 options: for the 3d-effect and an unbroken connection between the start/end scales
    frame $tw.fr_opt.f1 -borderwidth 1 -relief sunken
    frame $tw.fr_opt.f2
    frame $tw.fr_opt.f3
    frame $tw.fr_opt.f4

    # define two scales to choose the start and end time of the time vector
    # the values (variables) of the sliders equal the relativ x-positions of the sliders
    # MoveSliderStart/End: update the entry fields with the new start/end values
    # all options are necessary for an unbroken conncetion between the two scales
    # bind commands existing!
    scale $tw.fr_opt.f1.scale1 -from 10 -to 336 -orient horizontal -variable FC::Res(SliderStart) \
                               -borderwidth 0 -sliderlength 20 -length 346 \
                               -comm "FC::MoveSliderStart" -showvalue 0 -highlightthickness 0
    scale $tw.fr_opt.f1.scale2 -from 10 -to 10 -orient horizontal -variable FC::Res(SliderEnd) \
                               -borderwidth 0 -sliderlength 20 -length 20 \
                               -comm "FC::MoveSliderEnd" -showvalue 0 -highlightthickness 0

    gridrc 1 1 $tw.fr_opt.f1.scale1
    gridrc 1 2 $tw.fr_opt.f1.scale2

    # define two entries to enter the start and end time of the time vector manually
    # bind commands existing!
    label $tw.fr_opt.f2.start -text [l "Start \[s]" "Start \[s]"]
    label $tw.fr_opt.f2.stop -text [l "End \[s]" "Ende \[s]"]
    entry $tw.fr_opt.f2.e_start -textvariable FC::Res(valueStart) -width 7 -justify right
    entry $tw.fr_opt.f2.e_stop -textvariable FC::Res(valueEnd) -width 7 -justify right
    gridrc 1 1 $tw.fr_opt.f2.start -sticky w -padx 3 -pady {9 0}
    gridrc 1 2 $tw.fr_opt.f2.e_start -sticky w -padx {18 0} -pady {9 0}
    gridrc 1 3 $tw.fr_opt.f2.stop -sticky e -pady {9 0}
    gridrc 1 4 $tw.fr_opt.f2.e_stop -sticky e -padx 3 -pady {9 0}

    # define an offset entry/scale to choose an offset for the time vector
    # bind command existing!
    label $tw.fr_opt.f3.offs -text "Offset \[s]"
    entry $tw.fr_opt.f3.e_offs -textvariable FC::Res(valueOffset) -width 7 -justify right
    scale $tw.fr_opt.f3.scale3 -from "-100" -to "100" -orient horizontal \
                               -variable FC::Res(SliderOffset) -sliderlength 20 -length 242 \
                               -comm "FC::MoveSliderOffset" -showvalue 0  -resolution 0.01
    gridrc 1 1 $tw.fr_opt.f3.offs -sticky w -padx 3 -pady {3 0}
    gridrc 1 2 $tw.fr_opt.f3.e_offs -sticky w -padx {11 3} -pady {3 0}
    gridrc 1 3 $tw.fr_opt.f3.scale3 -sticky w -pady {2 0}

    # define an entry for the step size of the time vector (initial step size comes from Input File!)
    # and define a Reset button the set all scales and entries to default values
    label $tw.fr_opt.f4.step -text [l "Step size \[s]    " "Schrittweite \[s]"]
    entry $tw.fr_opt.f4.e_step -textvariable FC::Res(stepSize) -width 4 -justify right
    button $tw.fr_opt.f4.reset -text "Reset" -comm FC::ResOptionsReset
    gridrc 1 1 $tw.fr_opt.f4.step -sticky wn -padx {3 0} -pady {3 0}
    gridrc 1 2 $tw.fr_opt.f4.e_step -sticky wn -padx {4 0} -pady {3 0}
    gridrc 1 4 $tw.fr_opt.f4.reset -sticky es -padx 4 -pady {2 3}

    # draw the Time Base options
    gridrc 1 1 $tw.fr_opt.f2 -sticky we -padx 0
    gridrc 2 1 $tw.fr_opt.f1 -sticky we -padx 5 -pady {3 0}
    grid columnconfigure $tw.fr_opt.f2 3 -weight 1
    gridrc 3 1 $tw.fr_opt.f3 -sticky we -padx 0
    grid columnconfigure $tw.fr_opt.f3 3 -weight 1
    gridrc 4 1 $tw.fr_opt.f4 -sticky we -padx 0
    grid columnconfigure $tw.fr_opt.f4 3 -weight 1

    # define a frame for the OutputQuantities (channels) which can be selected (tw row 6)
    FrameWithTitle $tw.fr_var [l "OutputQuantities" "OutputQuantities"]
    frame $tw.fr_var.f1
    frame $tw.fr_var.f2

    # define a listbox (+scrollbar) including all channels (which can be found in the result file)
    # selected channels are blue highlighted, deselected are black
    # all channles are listed in Res(allQuants), selected channels in Res(selQuants)
    listbox $tw.fr_var.f1.lb_var -width 25 -height 10 -selectmode single \
                                 -listvariable  FC::Res(allQuants) \
                                 -yscrollcommand "$tw.fr_var.f1.lb_scroll set" \
                                 -selectbackground white  -activestyle none \
                                 -selectforeground black -foreground blue
    # set the color of the initial listbox entry to black (first element)
    $tw.fr_var.f1.lb_var itemconfigure 0 -foreground black
    scrollbar $tw.fr_var.f1.lb_scroll -orient vertical -command [list $tw.fr_var.f1.lb_var yview]
    gridrc 1 1 $tw.fr_var.f1.lb_var -sticky w -padx {3 0} -pady {9 3}
    gridrc 1 2 $tw.fr_var.f1.lb_scroll -sticky news -pady {10 4}

    # define two buttons to select/deselect all available channels
    # selection/deselection influences the foreground color of the listbox's elements 
    button $tw.fr_var.f2.all -text [l "Select all" "Alle auswählen"] -comm "FC::SelAllQuants"
    button $tw.fr_var.f2.del -text [l " Clear all" "Alle abwählen  "] -comm "FC::ClearAllQuants"
    label $tw.fr_var.f2.l_sel -textvariable FC::Res(labelQuants)
    gridrc 1 1 $tw.fr_var.f2.all -sticky wn -padx 9 -pady {9 0}
    gridrc 2 1 $tw.fr_var.f2.del -sticky wn -padx 9 -pady 3
    gridrc 3 1 $tw.fr_var.f2.l_sel -sticky wn -padx 9

    # draw the OutputQuantities options
    gridrc 5 1 $tw.fr_var.f1 -sticky w
    gridrc 5 2 $tw.fr_var.f2 -sticky wn 

    # draw the whole result option tw GUI
    gridrc 5 1 $tw.fr_opt -sticky we -padx 4 -pady 10 -columnspan 3
    gridrc 6 1 $tw.fr_var -sticky we -padx 4 -pady {0 5} -columnspan 3
    grid columnconfigure $tw.fr_var 3 -weight 1

### bindings for the result options
    # moving slider Start (start time) leads to a re-calculation of the options (of both scales!)
    bind $tw.fr_opt.f1.scale1 <B1-Motion> {
        # was slider End moved before? 
        if {$FC::Res(movedLast)=="SliderEnd"} {
            # options of Slider Start have to be adapted once, before moving it!
            FC::adaptScaleStart $FC::tw.fr_opt.f1.scale1 $FC::tw.fr_opt.f1.scale2
            break
        }
    }

    # adapt slider Start for manual entering the start time
    bind $tw.fr_opt.f2.e_start <ButtonRelease-1> {
        # was slider End moved before?
        if {$FC::Res(movedLast)=="SliderEnd"} {
            FC::adaptScaleStart $FC::tw.fr_opt.f1.scale1 $FC::tw.fr_opt.f1.scale2
        }
    }

    # move slider End for end time
    bind $tw.fr_opt.f1.scale2 <B1-Motion> {
        # was slider Start moved before?
        if {$FC::Res(movedLast)=="SliderStart"} {
            FC::adaptScaleEnd $FC::tw.fr_opt.f1.scale1 $FC::tw.fr_opt.f1.scale2
            break
        }
    }

    # adapt slider End for manual entering the end time
    bind $tw.fr_opt.f2.e_stop <ButtonRelease-1> {
        # was slider Start moved before?
        if {$FC::Res(movedLast)=="SliderStart"} {
            FC::adaptScaleEnd $FC::tw.fr_opt.f1.scale1 $FC::tw.fr_opt.f1.scale2
        }
    }

    # calculate the slider Start position from the start entry
    bind $tw.fr_opt.f2.e_start <Return> {
        # limit the range (to the minimum)
        if {$FC::Res(valueStart) < $FC::Res(limitStart)} {
           set FC::Res(valueStart) $FC::Res(limitStart)
        }
        # limit the range to slider End position
        if {$FC::Res(valueStart) > $FC::Res(valueEnd)} {
            set FC::Res(valueStart) $FC::Res(valueEnd)
        }
        # calculate the slider Start position (linear relation between slider and entry value)
        set FC::Res(SliderStart) [expr 328.0/($FC::Res(limitEnd)-$FC::Res(limitStart)) *\
             ($FC::Res(valueStart)-$FC::Res(limitStart))+10]
    }

    # calculate the slider End position from the end entry
    bind $tw.fr_opt.f2.e_stop <Return> {
        # limit the range (to the maximum)
        if {$FC::Res(valueEnd) > $FC::Res(limitEnd)} {
            set FC::Res(valueEnd) $FC::Res(limitEnd)
        }
        # limit the range to slider Start position
        if {$FC::Res(valueEnd) < $FC::Res(valueStart)} {
            set FC::Res(valueEnd) $FC::Res(valueStart)
        }
        # calculate the slider End position (linear relation between slider and entry value)
        set FC::Res(SliderEnd) [expr 326.0/($FC::Res(limitEnd)-$FC::Res(limitStart)) *\
             ($FC::Res(valueEnd)-$FC::Res(valueStart))+10]
    }

    # calculate the slider Offset position from the offset entry
    bind $tw.fr_opt.f3.e_offs <Return> {

        if {$FC::Res(valueOffset)>[$FC::tw.fr_opt.f3.scale3 cget -to]} {
            $FC::tw.fr_opt.f3.scale3 configure -to $FC::Res(valueOffset)
            $FC::tw.fr_opt.f3.scale3 configure -from -$FC::Res(valueOffset)
        }
        if {$FC::Res(valueOffset)<[$FC::tw.fr_opt.f3.scale3 cget -from]} {
            $FC::tw.fr_opt.f3.scale3 configure -from $FC::Res(valueOffset)
            $FC::tw.fr_opt.f3.scale3 configure -to [expr $FC::Res(valueOffset)*-1]
        }
        set FC::Res(SliderOffset) $FC::Res(valueOffset)
    }

    # select the OutputQuantities and adapt the listbox (grey/black highlighted)+label
    bind $tw.fr_var.f1.lb_var <ButtonRelease-1> {
        if {[file exists $FC::Fc(InFile)]} {
            # get the active element from the listbox
            set i [$FC::tw.fr_var.f1.lb_var curselection]
            # findout the index of the quantities which are already selected
            set tmp_index [lsearch $FC::Res(selQuants) [lindex $FC::Res(allQuants) $i]]
            # no index found --> select this quantity (mark it blue)
            if { $tmp_index == -1 } {
                lappend FC::Res(selQuants) [lindex $FC::Res(allQuants) $i]
                $FC::tw.fr_var.f1.lb_var itemconfigure $i -foreground blue
                $FC::tw.fr_var.f1.lb_var configure -selectforeground blue
            # index found --> deselect this quantity (mark it black)
            }  else {
                set FC::Res(selQuants) [lreplace $FC::Res(selQuants) $tmp_index $tmp_index]
                $FC::tw.fr_var.f1.lb_var itemconfigure $i -foreground black
                $FC::tw.fr_var.f1.lb_var configure -selectforeground black
            }
            # update the label of selected/all quantities
            FC::setLabelQuants
        }
    }
### end of bindings
}
### end of CreateResOptions


### open and handle the file browser to select the input result file
proc OpenResFile {} {
    variable Fc
    variable Res

    set types {{{IPG result files} .erg} {{ASCII files} .dat}}
    set tmp_input [tk_getOpenFile -filetypes $types -initialfile $Fc(InFile)\
                                  -initialdir "[pwd]/SimOutput"\
                                  -title [l "Select result file..." "Wähle Ergebnisdatei aus..."]]
    # no file selected --> do not change the input file name
    if {$tmp_input!=""} {
        set Fc(InFile) $tmp_input
    } else {
        return
    }

    if {$Fc(OutFile) == "" || [file extension $Fc(InFile)] == [file extension $Fc(OutFile)]} {
        # set the ouput file name (=input file name)
        set Fc(OutFile) [file rootname $Fc(InFile)]
        # define the output file extension: .erg-->.dat / .dat-->.erg
        if {[file extension $Fc(InFile)] == ".erg"} {
            append Fc(OutFile) ".dat"
        } elseif {[file extension $Fc(InFile)] == ".dat"} {
            append Fc(OutFile) ".erg"
        }
    }

    # reading the input file and save the quantities into a list
    if {[file exists $Fc(InFile)]} {
        set Res(allQuants) ""
        set infoFile [open "$Fc(InFile).info"]

        while {1} {
            set line [gets $infoFile]
            if [regexp {File\.At\.[0-9]*\.Name} $line] {
                lappend Res(allQuants) [lindex $line 2]
            }
            if {[eof $infoFile]} {
                close $infoFile
                break
            }
        }

        # create a converted tmp_file.dat to filter out the time base
        set cmd $Fc(pathRes)
        lappend cmd "-s"
        lappend cmd "Time"
        lappend cmd "-o"
        lappend cmd "tmp_file"
        lappend cmd "-om"
        lappend cmd "ascii"
        lappend cmd $Fc(InFile)

        if [catch [exec {*}$cmd] err] {
            # write error messages into window
            tk_messageBox -message $err
        }

        set tmp_file [open "tmp_file.dat"]
        set Res(timeBase) ""
        # save time base into list
        foreach line2 [split [read $tmp_file] \n] {
            lappend Res(timeBase) [lindex $line2 0]
        }
        close $tmp_file
        file delete tmp_file.dat
        file delete tmp_file.dat.info
        # draw new (reset) option GUI --> make file parameters visible
        ResOptionsReset
    }
    # scroll to end positions
    .fc.fio.in xview end
    .fc.fio.out xview end
}
### end of OpenResFile


### open and handle the file browser to select the output file
proc SaveResFile {} {
    variable Fc

    if {[file extension $Fc(InFile)] == ".erg"} {
        set types {{{ASCII Files} .dat}}
    } elseif {[file extension $Fc(InFile)] == ".dat"} {
        set types {{{IPG Result Files} .erg}}
    } else {
        set types {{{IPG result files} .erg} {{ASCII files} .dat}}
    }
    set tmp_output [tk_getSaveFile -filetypes $types -initialfile $Fc(OutFile)\
                                   -initialdir "[pwd]/SimOutput" \
                                   -title [l "Select storage path..." "Wähle Speicherpfad aus..."]]
    # no file selected --> do not change the output file name
    if {$tmp_output!=""} {
        set Fc(OutFile) $tmp_output
    } else {
        return
    }
    .fc.fio.out xview end
}


### re-configure slider Start (S1) and slider End (S2) to move slider Start
proc adaptScaleStart {sc_start sc_end} {
    variable Res

    # S1 length = S1 Xposition + S2 Xposition
    $sc_start configure -length [expr $Res(SliderStart) + $Res(SliderEnd)]
    # S2 length = constant (incl. window width, paddings...) - S1 length
    $sc_end configure -length [expr 366 - [$sc_start cget -length]]
    # "-to"-value = constant (incl. window width, paddings...) - S2 length
    $sc_start configure -to [expr 366 - [$sc_end cget -length] - 10]
    $sc_end configure -to [expr [$sc_end cget -length] - 10]
    # remember the current end value
    set Res(tmpSlider) $Res(valueEnd)
    # move slider End to the left end
    set Res(SliderEnd) 10
    set Res(movedLast) "SliderStart"
}


### re-configure slider Start (S1) and slider End (S2) to move slider End
proc adaptScaleEnd {sc_start sc_end} {
    variable Res

    # remember the current length of scale end
    set tmp [$sc_end cget -length]
    # calculate length of scale End depending on window width and slider Start position
    $sc_end configure -length [expr 366 - $Res(SliderStart) - 10]
    $sc_end configure -length [$sc_end cget -length]
    # calculate length of scale Start depending on slider Start position
    $sc_start configure -length [expr $Res(SliderStart) + 10]
    # calculate the "-to"-value (depending on slider Start position)
    $sc_start configure -to $Res(SliderStart)
    # calculate the "-to"-value (depending on scale End length)
    $sc_end configure -to [expr [$sc_end cget -length] - 10]
    # calculate the value of slider End depending on the length of scale End
    # and the remembered length of scale End
    set Res(SliderEnd) [expr [$sc_end cget -length] - $tmp + 10]
    ### set the remembered end value
    set Res(valueEnd) $Res(tmpSlider)
    $FC::tw.fr_opt.f2.e_stop configure -textvariable FC::Res(valueEnd)
    set Res(movedLast) "SliderEnd"
}


### calculate the value for the start entry from the slider Start position
proc MoveSliderStart {value} {
    variable Res

    # calculate the start value belonging to slider Start (linear relation between slider and entry)
    set Res(valueStart) [expr ($Res(limitEnd)-$Res(limitStart))/326.0*($Res(SliderStart)-10)+$Res(limitStart)]
    # round to 2 decimal places
    set Res(valueStart) [format %.2f [expr round($Res(valueStart)*100)/100.0]]
    set Res(valueEnd) $Res(tmpSlider)
}


### calculate the value for the end entry from the slider End position
proc MoveSliderEnd {value} {
    variable Res

    # calculate the end value belonging to slider End (linear relation between slider and entry)
    set Res(valueEnd) [expr ($Res(limitEnd)-$Res(limitStart))/326.0 * ($Res(SliderEnd)-10)+$Res(limitStart)+$Res(valueStart)-$Res(limitStart)]
    # round to 2 decimal places
    set Res(valueEnd) [format %.2f [expr round($Res(valueEnd)*100)/100.0]]

    # reload end time after slider 2 length adaption (when moving slider 1)
    if {$Res(movedLast)=="SliderStart"} {
        set Res(valueEnd) $Res(tmpSlider)
    }

    # set end time to end limit once (initial value)
    if {$Res(initialFlag)} {
        set Res(valueEnd) $Res(limitEnd)
        set Res(initialFlag) 0
    }
}


### calculate the value for the offset entry from the slider Offset position
proc MoveSliderOffset {value} {
    variable Res

    set Res(valueOffset) $Res(SliderOffset)
    # round to 2 decimal places
    set Res(valueOffset) [format %.2f [expr round($Res(valueOffset)*100)/100.0]]
}


### select all OutputQuantities in listbox
proc SelAllQuants {} {
    variable tw
    variable Fc
    variable Res

    if {[file exists $Fc(InFile)]} {
        for {set i 0} {$i < [$tw.fr_var.f1.lb_var size]} {incr i 1} {
            $tw.fr_var.f1.lb_var itemconfigure $i -foreground blue
            $tw.fr_var.f1.lb_var configure -selectforeground blue
        }
        set Res(selQuants) $Res(allQuants)
        setLabelQuants
    }
}


### deselect all OutputQuantities in listbox
proc ClearAllQuants {} {
    variable tw
    variable Fc
    variable Res

    if {[file exists $Fc(InFile)]} {
        for {set i 0} {$i < [$tw.fr_var.f1.lb_var size]} {incr i 1} {
            $tw.fr_var.f1.lb_var itemconfigure $i -foreground black
            $tw.fr_var.f1.lb_var configure -selectforeground black
        }
        set Res(selQuants) ""
        setLabelQuants
    }
}


### reset the result option GUI to the default values
proc ResOptionsReset {} {
    variable tw
    variable Fc
    variable Res

    # default settings for scales
    set Res(SliderStart) 10
    set Res(SliderEnd) 10
    set Res(movedLast) "SliderStart"
    set Res(initialFlag) 1
    set Res(SliderOffset) 0
    set Res(valueOffset) 0.0

    # input file loaded? 
    if {[file exists $Fc(InFile)]} {
        # set (update) the limits and values according to the loaded file
        # start time (first time element in vector)
        set Res(limitStart) [lindex $Res(timeBase)  3]
        # end time (last time element in vector)
        set Res(limitEnd) [lindex $Res(timeBase) [expr [llength $Res(timeBase)]-2] ]
        set Res(valueStart) $Res(limitStart)
        set Res(valueEnd) $Res(limitEnd)
        set Res(tmpSlider) $Res(limitEnd)
        # set (update) the step size, calculated (mean value) from the time vector
        # 4 time elements needed
        set Res(stepSizeFile) [expr ([lindex $Res(timeBase) 6]-[lindex $Res(timeBase) 5]+\
                              [lindex $Res(timeBase) 4]-[lindex $Res(timeBase) 3])/2]
        # round to 3 decimal places
        set Res(stepSizeFile) [expr round($Res(stepSizeFile)*1000)/1000.0]
        # the default step size for conversion is the step size from the file
        set Res(stepSize) $Res(stepSizeFile)
        # select all (available) quantities and update the label (selected/all)
        set Res(selQuants) $Res(allQuants)
        setLabelQuants
    } else {
        # set (update) the limits and values according to the default settings
        set Res(limitStart) 0.0
        set Res(limitEnd) 1000.00
        set Res(valueStart) "$Res(limitStart)"
        set Res(valueEnd) "$Res(limitEnd)"
        set Res(tmpSlider) "$Res(limitEnd)"
        set Res(labelQuants) [l "(0 from 0 selected)" "(0 von 0 ausgewählt)"]
        set Res(stepSize) 0.02
        set Res(selQuants) ""
        set Res(allQuants) [l "{Select the channels here!}" "{Selektiere hier die Kanäle!}"]
    }
    # result option GUI existent?
    if {[winfo exists $tw.fr_opt] && $Fc(options)!="hidd"} {
        # destroy result option GUI
        destroy $tw.fr_var
        destroy $tw.fr_opt
        destroy $tw.fr_var_title
        destroy $tw.fr_opt_title
        set Fc(options) "dest"
        # draw a new result option GUI after destroying...
        if {$Fc(Type)=="res"} {
            CreateAddOpt
        }
    }
}
### end of ResOptionsReset


### update/set the number of selected QutputQuantities
proc setLabelQuants {} {
    set FC::Res(labelQuants) [l "([llength $FC::Res(selQuants)] from [llength $FC::Res(allQuants)]\
                              selected)" "([llength $FC::Res(selQuants)] von \
                              [llength $FC::Res(allQuants)] ausgewählt)"]
}


### start the conversion with resutil
proc convertResult {} {
    variable Fc
    variable Res

        if {$Res(selQuants) == ""} {
            puts stderr "Error!\nPlease define OutputQuantities..."
            tk_messageBox -message "Error!\nPlease define OutputQuantities..."
            return
        }
        # define output name --> delete extension (.erg/.dat) from string
        # correct extension is automatically added
        set cmd "$Fc(pathRes)"
	lappend cmd "-o"
	lappend cmd "[file rootname $Fc(OutFile)]"
        # define output format
        if {[file extension $Fc(InFile)] == ".erg"} {
            lappend cmd "-om"
            lappend cmd "ascii"
        } elseif {[file extension $Fc(InFile)] == ".dat"} {
            lappend cmd "-om"
            lappend cmd "erg"
        }
        lappend cmd "-start"
        lappend cmd "$Res(valueStart)"
        lappend cmd "-end"
        lappend cmd "$Res(valueEnd)"
        # offset existent?
        if {$Res(valueOffset)!=0} {
            lappend cmd "-toff"
            lappend cmd "$Res(valueOffset)"
        }

        # was the step size (via entry) changed?
        if {$Res(stepSize) != $Res(stepSizeFile)} {
            lappend cmd "-dtmin" ;# add the entered step size
            lappend cmd "$Res(stepSize)"
        }
        set Channels ""
        # all quantities selected?
        if {[llength $Res(selQuants)]!=[llength $Res(allQuants)]} {
            for {set i 0} {$i < [llength $Res(selQuants)]} {incr i 1} {
                if {$i == [expr [llength $Res(selQuants)]-1]} {
                    # set the last output channel
                    append Channels "[lindex $Res(selQuants) $i]"
                } else {
                    # set the output channels (comma seperated)
                    append Channels "[lindex $Res(selQuants) $i],"
                }
            }
            # add the channels to the string
            lappend cmd "-s"
            lappend cmd "$Channels"
        } else {
            # or add all channels
            lappend cmd "-a"
        }
        # define the input file name
        lappend cmd "$Fc(InFile)"

    # execute resutil
    if [catch [exec {*}$cmd] err] {
        # write error messages into window
        tk_messageBox -message $err

    } else {
        tk_messageBox -message [l "File created successfully!" "Datei erfolgreich konvertiert!"]
    }
}
### end of convertResult

}
### end of namespace FC

