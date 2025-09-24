foreach name  {NaCl_Box_50 CA_15mM_1264_Box_50 MG_15mM_1264_Box_50} {
mol load parm7 ../${name}_rw.prmtop 
mol addfile ../${name}_run1_rw.nc first 2500 last 12500 step 25 waitfor all
mol addfile ../${name}_run2_rw.nc first 2500 last 12500 step 25 waitfor all
mol addfile ../${name}_run3_rw.nc first 2500 last 12500 step 25 waitfor all

source rename_chainID_with_tail.tcl
set nframes [expr  [molinfo top get numframes] - 1 ]

set outfile0 [open ./dna_tail/dna_tail_mean_contacts_${name}_all.dat w]
set outfile1 [open ./dna_tail/dna_tail_mean_contacts_${name}_h3.dat w]
set outfile2 [open ./dna_tail/dna_tail_mean_contacts_${name}_h4.dat w]
set outfile3 [open ./dna_tail/dna_tail_mean_contacts_${name}_h2a.dat w]
set outfile4 [open ./dna_tail/dna_tail_mean_contacts_${name}_h2b.dat w]
set outfile5 [open ./dna_tail/dna_tail_mean_contacts_${name}_h2a_N.dat w]
set outfile6 [open ./dna_tail/dna_tail_mean_contacts_${name}_h2a_C.dat w]
	
puts -nonewline $outfile0 "Time"
puts -nonewline $outfile1 "Time"
puts -nonewline $outfile2 "Time"
puts -nonewline $outfile3 "Time"
puts -nonewline $outfile4 "Time"
puts -nonewline $outfile5 "Time"
puts -nonewline $outfile6 "Time"

for { set r1 -73 } { $r1<=73 } { incr r1 } {
puts -nonewline $outfile0 "\t$r1"
puts -nonewline $outfile1 "\t$r1"
puts -nonewline $outfile2 "\t$r1"
puts -nonewline $outfile3 "\t$r1"
puts -nonewline $outfile4 "\t$r1"
puts -nonewline $outfile5 "\t$r1"
puts -nonewline $outfile6 "\t$r1"
}
puts -nonewline $outfile0 "\n"
puts -nonewline $outfile1 "\n"
puts -nonewline $outfile2 "\n"
puts -nonewline $outfile3 "\n"
puts -nonewline $outfile4 "\n"
puts -nonewline $outfile5 "\n"
puts -nonewline $outfile6 "\n"
for { set i 0 } { $i<=$nframes } { incr i } {
set time [expr 1 * $i]
puts -nonewline $outfile0 [format "%.3f" "$time"]
puts -nonewline $outfile1 [format "%.3f" "$time"]
puts -nonewline $outfile2 [format "%.3f" "$time"]
puts -nonewline $outfile3 [format "%.3f" "$time"]
puts -nonewline $outfile4 [format "%.3f" "$time"]
puts -nonewline $outfile5 [format "%.3f" "$time"]
puts -nonewline $outfile6 [format "%.3f" "$time"]
puts  [format "Time %.3f" "$time"]
unset time


for { set r1 -73 } { $r1<=73 } { incr r1 } {
set sel0 [atomselect top "((segname CHA CHE and resid 1 to 37) or (segname CHB  CHF and resid 1 to 20) or (segname CHC CHG and resid 1 to 13 118 to 129) or (segname CHD CHH and resid 1 to 26)) and noh"  frame $i]
set sel1 [atomselect top "(segname CHA CHE and resid 1 to 37) and noh"  frame $i]
set sel2 [atomselect top "(segname CHB CHF and resid 1 to 20) and noh"  frame $i]
set sel3 [atomselect top "(segname CHC CHG and resid 1 to 13 118 to 129) and noh"  frame $i]
set sel4 [atomselect top "(segname CHD CHH and resid 1 to 26) and noh"  frame $i]
set sel5 [atomselect top "(segname CHC CHG and resid 1 to 13 ) and noh"  frame $i]
set sel6 [atomselect top "(segname CHC CHG and resid 118 to 129) and noh"  frame $i]

set r2 [expr $r1 * (-1)]
set sel7 [atomselect top "((segname CHI and resid '$r1') or (segname CHJ and resid '$r2')) and noh" frame $i]
set contacts0 [lindex [measure contacts  4.0 $sel0 $sel7] 1]
set contacts1 [lindex [measure contacts  4.0 $sel1 $sel7] 1]
set contacts2 [lindex [measure contacts  4.0 $sel2 $sel7] 1]
set contacts3 [lindex [measure contacts  4.0 $sel3 $sel7] 1]
set contacts4 [lindex [measure contacts  4.0 $sel4 $sel7] 1]
set contacts5 [lindex [measure contacts  4.0 $sel5 $sel7] 1]
set contacts6 [lindex [measure contacts  4.0 $sel6 $sel7] 1]

set nc0 [llength $contacts0]
set nc1 [llength $contacts1]
set nc2 [llength $contacts2]
set nc3 [llength $contacts3]
set nc4 [llength $contacts4]
set nc5 [llength $contacts5]
set nc6 [llength $contacts6]

puts -nonewline $outfile0 "\t$nc0"
puts -nonewline $outfile1 "\t$nc1"
puts -nonewline $outfile2 "\t$nc2"
puts -nonewline $outfile3 "\t$nc3"
puts -nonewline $outfile4 "\t$nc4"
puts -nonewline $outfile5 "\t$nc5"
puts -nonewline $outfile6 "\t$nc6"

$sel0 delete
$sel1 delete
$sel2 delete
$sel3 delete
$sel4 delete
$sel5 delete
$sel6 delete
$sel7 delete

unset contacts0 
unset contacts1 
unset contacts2 
unset contacts3 
unset contacts4 
unset contacts5 
unset contacts6

unset nc0 
unset nc1 
unset nc2 
unset nc3
unset nc4
unset nc5 
unset nc6 

}
puts -nonewline $outfile0 "\n"
puts -nonewline $outfile1 "\n"
puts -nonewline $outfile2 "\n"
puts -nonewline $outfile3 "\n"
puts -nonewline $outfile4 "\n"
puts -nonewline $outfile5 "\n"
puts -nonewline $outfile6 "\n"
}

close $outfile0
close $outfile1 
close $outfile2
close $outfile3
close $outfile4
close $outfile5
close $outfile6
}
exit

