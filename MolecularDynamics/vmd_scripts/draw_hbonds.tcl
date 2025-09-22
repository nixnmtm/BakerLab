# ---------- draw_hbonds.tcl ----------
# You set: molid, HBONDS list, then call draw_hbonds_from_list
# Each HBONDS item: {d_resn d_resid d_atom  a_resn a_resid a_atom  occupancy}

# Example color bands and radii
set rmin 0.10      ;# radius at 0 occupancy (Å in VMD screen units)
set rmax 0.60      ;# radius at 100% occupancy
set lowThr  0.20   ;# 0–20% = low
set midThr  0.60   ;# 20–60% = mid; >60% = high
set colLow  0      ;# VMD color ID (0=blue)
set colMid  3      ;# 3=orange
set colHigh 1      ;# 1=red
set labelOcc 1     ;# 1 to draw tiny occupancy labels near the bond

proc _hb_color_for_occ {occ lowThr midThr colLow colMid colHigh} {
    if {$occ <= $lowThr} { return $colLow }
    if {$occ <= $midThr} { return $colMid }
    return $colHigh
}

proc draw_hbonds_from_list {molid HBONDS rmin rmax lowThr midThr colLow colMid colHigh labelOcc} {
    draw color white
    foreach item $HBONDS {
        # unpack: d_resn d_resid d_atom  a_resn a_resid a_atom  occupancy
        lassign $item dresn dresi datm aresn aresi aatm occ
        # selections (single atoms)
        set sel1 [atomselect $molid "resname $dresn and resid $dresi and name $datm"]
        set sel2 [atomselect $molid "resname $aresn and resid $aresi and name $aatm"]
        if {[$sel1 num] < 1 || [$sel2 num] < 1} {
            $sel1 delete; $sel2 delete
            continue
        }
        set p1 [lindex [$sel1 get {x y z}] 0]
        set p2 [lindex [$sel2 get {x y z}] 0]
        $sel1 delete
        $sel2 delete

        # radius scaled by occupancy (0..1)
        set r [expr {$rmin + ($rmax - $rmin) * $occ}]
        set col [_hb_color_for_occ $occ $lowThr $midThr $colLow $colMid $colHigh]

        draw color $col
        draw cylinder $p1 $p2 radius $r filled yes resolution 12

        if {$labelOcc} {
            # place a tiny label at 40% along the cylinder
            set x [expr {[lindex $p1 0]*0.6 + [lindex $p2 0]*0.4}]
            set y [expr {[lindex $p1 1]*0.6 + [lindex $p2 1]*0.4}]
            set z [expr {[lindex $p1 2]*0.6 + [lindex $p2 2]*0.4}]
            draw color white
            draw text "$x $y $z" [format "%.0f%%" [expr {$occ*100.0}]] size 0.7
        }
    }
}

# ------------- HOW TO USE -------------
# 1) In VMD Tk Console:
#    source draw_hbonds.tcl
# 2) Set your molecule id:
#    set molid 0
# 3) Paste your pairs block here (from Python emit_vmd_pairs):
#    set HBONDS {
#      { LYS 184 {NZ}  ATP 1 {O3G}  0.92 }
#      { ASN 243 {HD22}  ATP 1 {O3B}  0.47 }
#      ...
#    }
# 4) Draw:
#    draw delete all
#    draw_hbonds_from_list $molid $HBONDS $rmin $rmax $lowThr $midThr $colLow $colMid $colHigh $labelOcc

