set system_basename peptoid
set replicate_ids [split $argv ","]
set sidelen 20. ;# half the box side length
set outprefix pad

foreach rep $replicate_ids {
    mol new ${system_basename}_rep${rep}_min.pdb
    # align S-S to Z axis
    set all [atomselect top all]
    set ss_pos [[atomselect top "name S1"] get {x y z}]
    set ss_ax [vecnorm [vecsub [lindex $ss_pos 1] [lindex $ss_pos 0]]]
    $all move [transvecinv $ss_ax]
    $all move [transaxis y -90]
    $all writepdb ${system_basename}_rep${rep}_min_align.pdb
    $all delete

    # get padding dims
    mol new ${system_basename}_rep${rep}_min_align.pdb
    set all [atomselect top all]
    set lims [measure minmax $all]
    set maxs [vecsub "$sidelen $sidelen $sidelen" [lindex $lims 1]]
    set mins [vecscale -1. [vecsub "-$sidelen -$sidelen -$sidelen" [lindex $lims 0]]]
    set pad [vecscale 0.5 [vecadd $mins $maxs]]

    # save pad dims
    set fp [open ${outprefix}_rep${rep}.dat "w"]
    puts $fp $pad
    close $fp
}

exit
