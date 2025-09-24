proc renumber { sel start } {
	if { [$sel num] == 0 } {
		puts "Error in renumber: empty selection!"
		return
	}
	set oresid [ $sel get resid ]
	set delta [ expr $start - [ lindex $oresid 0] ]
	set nresid { }
	foreach r $oresid {
		lappend nresid [ expr $r + $delta ]
	}
	$sel set resid $nresid
}


set sel [atomselect top "resid 1 to 147 "]
$sel set chain I
$sel set segname CHI
renumber $sel -73
$sel delete

set sel [atomselect top "resid 148 to 294 "]
$sel set chain J
$sel set segname CHJ
renumber $sel -73
$sel delete

set sel [atomselect top "resid 295 to 423 "]
$sel set chain C
$sel set segname CHC
renumber $sel 1
$sel delete

set sel [atomselect top "resid 424 to 552 "]
$sel set chain G
$sel set segname CHG
renumber $sel 1
$sel delete


set sel [atomselect top "resid 553 to 677 "]
$sel set chain D
$sel set segname CHD
renumber $sel 1
$sel delete

set sel [atomselect top "resid 678 to 802 "]
$sel set chain H
$sel set segname CHH
renumber $sel 1
$sel delete

set sel [atomselect top "resid 803 to 937 "]
$sel set chain A
$sel set segname CHA
renumber $sel 1
$sel delete


set sel [atomselect top "resid 938 to 1072 "]
$sel set chain E
$sel set segname CHE
renumber $sel 1
$sel delete


set sel [atomselect top "resid 1073 to 1174 "]
$sel set chain B
$sel set segname CHB
renumber $sel 1
$sel delete


set sel [atomselect top "resid 1175 to 1276 "]
$sel set chain F
$sel set segname CHF
renumber $sel 1
$sel delete

