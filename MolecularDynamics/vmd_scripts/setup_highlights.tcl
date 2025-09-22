# 00_setup_highlights.tcl
# Usage (after loading files):
#   source 00_setup_highlights.tcl
#   setup_scene top WT      ;# or: setup_scene top R108W

proc setup_scene {molid variant} {
    # --- display tuning for smooth on-screen + consistent renders ---
    display rendermode GLSL
    display ambientocclusion on
    display shadows on
    display antialias on
    display depthcue off
    axes location off
    color Display Background white
    # higher AA for on-screen (export uses renderer AA too)
    display aaquality 8

    # --- selections (SM-B2 numbering) ---
    # Shared (WT canonical contacts)
    set sel_atp  [atomselect $molid "resname ATP"]
    set sel_mg   [atomselect $molid "name MG"]
    set sel_prot [atomselect $molid "protein"]

    # WT “canonical” cast
    set wt_resid_list {126 129 180 181 182 183 184 185 186 243 245 246 247}
    # R108W “rewired” cast (keep WT + Switch II glycine 469)
    set mut_resid_list {126 129 180 181 182 183 184 185 186 243 245 246 247 469}

    if {$variant eq "WT"} {
        set keyresids $wt_resid_list
    } else {
        set keyresids $mut_resid_list
    }
    set sel_keys  [atomselect $molid "protein and resid [join $keyresids " "]"]

    # --- clear reps and build consistent stack ---
    mol top $molid
    set nrep [molinfo $molid get numreps]
    for {set i 0} {$i < $nrep} {incr i} { mol delrep 0 $molid }

    # 1) protein cartoon (faded)
    mol representation NewCartoon
    mol selection "protein"
    mol color ColorID 8
    mol addrep $molid

    # 2) ATP sticks
    mol representation Licorice 0.25 12 12
    mol selection "resname ATP"
    mol color Name
    mol addrep $molid

    # 3) Mg2+ big sphere
    mol representation VDW 2.2
    mol selection "name MG"
    mol color ColorID 2
    mol addrep $molid

    # 4) Key residues sticks (highlighted)
    mol representation Licorice 0.28 12 12
    mol selection "protein and resid [join $keyresids " "]"
    mol color ColorID 1
    mol addrep $molid

    # Optional: color WT vs R108W differently
    if {$variant eq "R108W"} {
        # give mutant key residues a distinct color
        set repidx [expr {[molinfo $molid get numreps] - 1}]
        mol modcolor $repidx $molid ColorID 4
    }

    # --- trajectory smoothing on visual reps (prettier motion) ---
    # increase smoothing window on all reps
    set reps [molinfo $molid get numreps]
    for {set r 0} {$r < $reps} {incr r} {
        # smoothrep: moving-average over N frames; 5–10 looks nice without lag
        mol smoothrep $molid $r 7
    }

    # --- optional: distance labels for marquee contacts per variant ---
    # choose a few to show; we’ll label one or two distances to reduce clutter
    # LYS184:NZ ↔ ATP:O3G (dominant in mutant; canonical in WT to β/γ)
    set sel_lys184 [atomselect $molid "protein and resid 184 and name NZ"]
    set sel_o3g    [atomselect $molid "resname ATP and name O3G"]
    if {[$sel_lys184 num] && [$sel_o3g num]} {
        label add Bonds [$sel_lys184 index] [$sel_o3g index]
    }

    # GLY469:N ↔ ATP:O1G/3G (mutant)
    if {$variant eq "R108W"} {
        foreach oname {O1G O3G} {
            set sel_g469n [atomselect $molid "protein and resid 469 and name N"]
            set sel_op    [atomselect $molid "resname ATP and name $oname"]
            if {[$sel_g469n num] && [$sel_op num]} {
                label add Bonds [$sel_g469n index] [$sel_op index]
            }
        }
    }

    # Asn126/ Lys129 adenine contacts (fade in WT; will appear long when broken)
    if {$variant eq "WT"} {
        set sel_asn126 [atomselect $molid "protein and resid 126 and name ND2"]
        set sel_n7     [atomselect $molid "resname ATP and name N7"]
        if {[$sel_asn126 num] && [$sel_n7 num]} {
            label add Bonds [$sel_asn126 index] [$sel_n7 index]
        }
        set sel_k129o  [atomselect $molid "protein and resid 129 and name O"]
        set sel_n6     [atomselect $molid "resname ATP and name N6"]
        if {[$sel_k129o num] && [$sel_n6 num]} {
            label add Bonds [$sel_k129o index] [$sel_n6 index]
        }
    }

    # Tidy labels
    label textsize 1.2
    label color Bonds black
    label size Bonds 0.8

    # Camera defaults for consistency
    display projection Orthographic
    rotate x by -10
    rotate y by 15
    scale by 1.2
    # lock to ATP centroid so the pocket stays centered (optional)
    set c [$sel_atp center]
    translate to $c
    # done
    puts "Scene setup complete for variant=$variant (molid=$molid)."
}

