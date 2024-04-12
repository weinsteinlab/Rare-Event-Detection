set readstartf [lindex $argv 0]
set outfile [lindex $argv 1]
set stride [lindex $argv 2]
set trajlength 1000

if { $readstartf > 0 } {
set delete_until [expr $readstartf - 1 ]
animate delete end $delete_until}
set nf [molinfo top get numframes]
echo "number of frames is $nf"

animate delete beg $trajlength
set nf [molinfo top get numframes]
echo "number of frames is $nf"

if {$stride > 1} {animate delete skip $stride}

set nf [molinfo top get numframes]
echo "number of frames is $nf"

set prot [atomselect top "protein"]
set rezz_list [$prot get resid]
set rezz_list [lsort -unique -integer $rezz_list]

set ff [open $outfile w]

for {set fr 0} {$fr < $nf} {incr fr} {
    set results []
    foreach rez $rezz_list {
        echo $rez
        set sel [atomselect top "protein and same residue as within 3.5 of (protein and resid $rez)" frame $fr]
        set contacts [$sel get resid]
        $sel delete
        set contacts [lsort -unique -integer $contacts]
        set results [concat $results "$contacts,"]
    }
    puts $ff $results
    unset results
    unset contacts
}


close $ff
exit
