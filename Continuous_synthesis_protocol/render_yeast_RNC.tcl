set traj_idx [lindex $argv 0]
set stage_idx [lindex $argv 1]
set outname [lindex $argv 2]

set traj_name "traj_"
append traj_name $traj_idx "_" $stage_idx ".pdb"

mol new $traj_name type pdb first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all

display resize 1200 900
display projection Orthographic
display rendermode GLSL
axes location off

color Display {Background} white

material add YangTrans1
material change ambient YangTrans1 0.000000
material change diffuse YangTrans1 0.730000
material change specular YangTrans1 0.600000
material change shininess YangTrans1 0.680000
material change mirror YangTrans1 0.000000
material change opacity YangTrans1 0.430000
material change outline YangTrans1 0.660000
material change outlinewidth YangTrans1 0.230000
material change transmode YangTrans1 1.000000

color Structure {Alpha Helix} purple
color Structure 3_10_Helix purple
color Structure Pi_Helix purple
color Structure Extended_Beta yellow
color Structure Bridge_Beta yellow
color Structure Turn cyan
color Structure Coil cyan

set sel [atomselect top "segname A"]
set chain_length [lindex [$sel get resid] end]

mol delrep 0 top
if { $stage_idx == 1 } {   
  if { $chain_length > 1 } {
    if { $chain_length > 10 } {
      mol representation NewCartoon 0.300000 24.000000 4.100000 0
    } else {
      mol representation Tube
    }
    #mol color ColorID 27
    mol color Structure
    set sel "segname A and resid 1 to "
    append sel [expr $chain_length - 1]
    mol selection $sel
    mol material AOEdgy
    mol addrep top
  }
  mol representation VDW
  mol color ColorID 27
  #mol color Structure
  set sel "segname A and resid "
  append sel $chain_length " and name CA"
  mol selection $sel
  mol material AOEdgy
  mol addrep top
  
} else {
  if { $chain_length > 10 } {
    mol representation NewCartoon 0.300000 24.000000 4.100000 0
  } else {
    mol representation Tube
  }
  #mol color ColorID 27
  mol color Structure
  mol selection "segname A"
  mol material AOEdgy
  mol addrep top
}

if { $stage_idx != 4 } { 
  mol representation QuickSurf 1.600000 1.000000 1.000000 3.000000
  mol color ColorID 3
  if { $stage_idx == 3 } { 
    mol selection "segname PtR"
  } else {
    mol selection "segname AtR PtR"
  }
  mol material EdgyShiny
  mol addrep top
}

if { $stage_idx == 1 } { 
  mol representation DynamicBonds 5.000000 1.000000 12.000000
  mol color ColorID 3
  set sel "segname AtR and resid 76 and name R or segname A and resid "
  append sel $chain_length " and name CA"
  mol selection $sel
  mol material EdgyShiny
  mol addrep top
  
  if { $chain_length > 1 } {
    mol representation DynamicBonds 5.000000 1.000000 12.000000
    mol color ColorID 3
    set sel "segname PtR and resid 76 and name R or segname A and resid "
    append sel [expr $chain_length - 1] " and name CA"
    mol selection $sel
    mol material EdgyShiny
    mol addrep top
  }
} elseif {$stage_idx == 2} {
  mol representation DynamicBonds 5.000000 1.000000 12.000000
  mol color ColorID 3
  set sel "segname AtR and resid 76 and name R or segname A and resid "
  append sel $chain_length " and name CA"
  mol selection $sel
  mol material EdgyShiny
  mol addrep top
} elseif {$stage_idx == 3} {
  mol representation DynamicBonds 5.000000 1.000000 12.000000
  mol color ColorID 3
  set sel "segname PtR and resid 76 and name R or segname A and resid "
  append sel $chain_length " and name CA"
  mol selection $sel
  mol material EdgyShiny
  mol addrep top
}

if { $stage_idx != 4 } { 
  mol representation QuickSurf 3.000000 2.000000 1.000000 3.000000
  mol color ColorID 8
  mol selection "not segname A AtR PtR and x<58"
  mol material YangTrans1
  mol addrep top

  mol representation QuickSurf 3.000000 2.000000 1.000000 3.000000
  mol color ColorID 8
  mol selection "not segname A AtR PtR and (x>=58 and y<20 and y>-20 and z>-15)"
  mol material YangTrans1
  mol addrep top

  mol representation QuickSurf 3.000000 2.000000 1.000000 3.000000
  mol color ColorID 8
  mol selection "(not segname A AtR PtR) and x>=58 and not (y<20 and y>-20 and z>-15)"
  mol material BrushedMetal
  #mol material YangTrans1
  mol addrep top
}

mol new /storage/home/yuj179/mygroup/ribosome/Yeast/add_tRNA/6q8y_cg/60S_tRNA_cg.psf type psf first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol addfile /storage/home/yuj179/mygroup/ribosome/Yeast/add_tRNA/6q8y_cg/60S_tRNA_cg.cor type cor first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all
mol delrep 0 top
mol representation QuickSurf 3.000000 2.000000 1.000000 3.000000
if { $stage_idx != 4 } {
  mol color ColorID 15
  mol selection {not (x<58 and x>-10 and y<32 and y > -32 and z>-5) and not x > 58}
  mol material Ghost
} else {
  mol color ColorID 8
  mol selection "all"
  mol material BrushedMetal
}
mol addrep top

#mol top 0
#display resetview
translate to -1 0 0
rotate y by -30 
scale by 1.2

render Tachyon $outname "$env(VMDDIR)/tachyon_LINUXAMD64 -aasamples 12 %s -format TARGA -trans_max_surfaces 1 -o %s.tga"

quit
