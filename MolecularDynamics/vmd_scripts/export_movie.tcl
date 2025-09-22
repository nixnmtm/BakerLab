# ============================================
# Robust movie exporter (WT/R108W), mac-friendly
# Usage examples:
#   export_movie top out_wt 54000 60000 20 60 snapshot
#   export_movie top out_mut 54000 60000 20 60 TachyonInternal
#   export_movie top out_wt 54000 60000 20 60   ;# (auto-pick renderer)
# ============================================
proc export_movie {molid outprefix start end step fps {renderer ""}} {
    # Normalize molid
    if {$molid eq "top"} { set molid [molinfo top get id] }
    if {$molid eq ""} { error "No molecule loaded / invalid molid" }

    # Validate frame range
    set nf [molinfo $molid get numframes]
    if {$nf <= 0} { error "Molecule $molid has no frames" }
    if {$end eq "last"} { set end [expr {$nf - 1}] }
    if {$start < 0} { set start 0 }
    if {$end >= $nf} { set end [expr {$nf - 1}] }
    if {$step < 1} { set step 1 }

    # Pick renderer if not provided
    set rlist [render list]
    if {$renderer eq ""} {
        if {[lsearch -exact $rlist "snapshot"] >= 0} {
            set renderer snapshot
        } elseif {[lsearch -exact $rlist "TachyonInternal"] >= 0} {
            set renderer TachyonInternal
        } else {
            set renderer [lindex $rlist 0]
        }
    }
    puts "export_movie: molid=$molid frames=$start..$end step=$step fps=$fps renderer=$renderer"
    puts "export_movie: available renderers: $rlist"

    # Simple, fast display defaults
    display projection Orthographic
    display resize 1920 1080
    display aaquality 4
    display ambientocclusion off
    display shadows off
    display depthcue off
    color Display Background white

    # Try to sanitize any malformed rep selections so render won't bail
    set nrep [molinfo $molid get numreps]
    for {set r 0} {$r < $nrep} {incr r} {
        if {[catch {mol repselect $molid $r} sel]} { continue }
        if {[catch {mol modselect $r $molid $sel} err]} {
            puts "WARN: rep $r had invalid selection '$sel' → using 'all' temporarily"
            mol modselect $r $molid "all"
        }
    }

    # Output directory
    set outdir "${outprefix}_frames"
    file mkdir $outdir
    set pad 5

    # Render loop
    set count 0
    for {set f $start} {$f <= $end} {incr f $step} {
        animate goto $f
        set fn [format "%s/%s_%0${pad}d.tga" $outdir $outprefix $count]
        if {[catch {render $renderer $fn} rerr]} {
            puts "WARN: render failed at frame $f with '$rerr' — continuing"
            continue
        }
        incr count
    }
    puts "export_movie: wrote $count frames → $outdir"
    puts "ffmpeg (mac HW encode): ffmpeg -y -framerate $fps -i ${outdir}/${outprefix}_%0${pad}d.tga -pix_fmt yuv420p -c:v h264_videotoolbox -b:v 12M ${outprefix}.mp4"
}

