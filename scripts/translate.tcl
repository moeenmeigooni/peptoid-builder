# 2024 Jul 23 Moeen Meigooni
#
package require pbctools

set system_basename peptoid
set replicate_ids [split $argv ","]

foreach rep $replicate_ids {
    # load 
    mol new ${system_basename}_water_rep${rep}.pdb

    # translate coords for OpenMM-style pbc
    set all [atomselect top all]
    set mm [measure minmax $all]
    set newbox [vecsub [lindex $mm 1] [lindex $mm 0]]
    $all moveby [vecscale -1.0 [lindex $mm 0]]

    # set pbc
    pbc set [list [lindex $newbox 0] [lindex $newbox 1] [lindex $newbox 2] 90. 90. 90.]

    # save inplace
    $all writepdb ${system_basename}_water_rep${rep}.pdb
}

exit
