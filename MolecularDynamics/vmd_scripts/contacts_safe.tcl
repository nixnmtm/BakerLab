# wt_atp_contacts_safe.tcl
# Show ONLY residues that H-bond (<3.5 Å, 30°) to ATP while you step frames.
# Also includes an optional frame loop to make frames for a short movie.

namespace eval hbvid {
    variable molid
    variable keyres {126 129 180 181 182 183 184 185 186 243 245 246 247 469}
    variable cutoff 3.5
    variable angle 30
    variable rep_dynamic -1
    variable sel_show
    variable sel_atp
    variable sel_prot
}

proc hbvid::setup {molid_in} {
    variable molid; set molid $molid_in
    variable sel_atp;  variable sel_prot; variable sel_show; variable rep_dynamic

    # Selections
    set sel_prot [atomselect $molid "protein"]
    set sel_atp  [atomselect $molid "resname ATP"]
    if {[$sel_atp num] == 0} {
        puts "WARN: No atoms for resname ATP. If ligand is ANP/ATP3, edit the selection."
    }

    # Clean reps
    mol top $molid
    set nrep [molinfo $molid get numreps]
    for {set i 0} {$i < $nrep} {incr i} { mol delrep 0 $molid }

    # 1) protein faint cartoon
    mol representation NewCartoon
    mol selection "protein"
    mol color ColorID 7
    mol material Glass1
    mol addrep $molid

    # 2) ATP sticks
    mol representation Licorice 0.25 12 12
    mol selection "resname ATP"
    mol color Name
    mol addrep $molid

    # 3) Mg2+ sphere
    mol representation VDW 0.5
    mol selection "name MG"
    mol color Name
    mol addrep $molid

    # 4) dynamic residue sticks (start empty)
    set sel_show [atomselect $molid "index -1"]
    mol representation Licorice 0.28 12 12
    mol selection "index -1"
    mol color ColorID 1
    mol addrep $molid
    set rep_dynamic [expr {[molinfo $molid get numreps] - 1}]

    # smoothing
    set reps [molinfo $molid get numreps]
    for {set r 0} {$r < $reps} {incr r} { mol smoothrep $molid $r 7 }

    puts "hbvid: setup complete on molid=$molid"
}

proc hb_update {} {
    global hb_molid hb_keyres hb_cutoff hb_angle hb_rep_dynamic hb_sel_show
    if {$hb_molid < 0} { puts "hb_update: run hb_setup first"; return }

    graphics $hb_molid delete all

    set protSel [atomselect $hb_molid "protein and resid [join $hb_keyres " "]"]
    set atpSel  [atomselect $hb_molid "resname ATP"]
    if {[$protSel num] == 0 || [$atpSel num] == 0} {
        $hb_sel_show set selection "index -1"
        mol modselect $hb_rep_dynamic $hb_molid "index -1"
        return
    }

    lassign [measure hbonds $hb_cutoff $hb_angle $protSel $atpSel] donIdx accIdx pairs

    # protein atoms to show
    set protContactIdx [lsort -unique $donIdx]
    if {[llength $protContactIdx] == 0} {
        $hb_sel_show set selection "index -1"
        mol modselect $hb_rep_dynamic $hb_molid "index -1"
        return
    }

    $hb_sel_show set selection "index [join $protContactIdx " "]"
    $hb_sel_show update
    mol modselect $hb_rep_dynamic $hb_molid "[$hb_sel_show text]"

    # draw bonds
    set n [llength $donIdx]
    for {set k 0} {$k < $n} {incr k} {
        set i [lindex $donIdx $k]
        set j [lindex $accIdx $k]
        set si [atomselect $hb_molid "index $i"]
        set sj [atomselect $hb_molid "index $j"]
        if {[$si num] && [$sj num]} {
            set ci [lindex [$si get {x y z}] 0]
            set cj [lindex [$sj get {x y z}] 0]
            graphics $hb_molid color 1
            graphics $hb_molid cylinder $ci $cj radius 0.08 filled yes
        }
        $si delete
        $sj delete
    }
    $protSel delete
    $atpSel delete
}


proc hbvid::rep_dynamic_id {} {
    # helper to fetch the (last) rep id reliably
    variable molid
    return [expr {[molinfo $molid get numreps] - 1}]
}

# Optional: frame loop to preview or render a SHORT segment (no Tk needed)
# Example: hbvid::playloop [molinfo top get id] 54000 60000 20
proc hbvid::playloop {molid start end step} {
    variable molid
    for {set f $start} {$f <= $end} {incr f $step} {
        animate goto $f
        hbvid::update
        display update
    }
}

