#!/usr/bin/perl

@l=(17,20,22,23,24,25,26,27,28,29,30,31,32,33,36,37,39,41,43,46);

$c = '/gpfs/home/epo2/work/software/charmm/executables/c35b5_deby_doubang_pull/charmm';
foreach $ele (@l)
  {
  $offset = $ele - 17 + 49;
  $cmd = "$c out=../shared_files/rnc_files/rnc_2adr_l"."$ele"."_ca psf=../shared_files/nascent_chains/2adr_l"."$ele"."_ca.psf crd=../shared_files/nascent_chains/2adr_l"."$ele"."_ca.cor prm=../shared_files/params/rnc-2adr_l"."$ele"."_nscal1.7_fnn1_go_bt.prm top=../shared_files/tops/rnc-2adr_l"."$ele"."_ca.top off=$offset < combine_2adr_and_ribosome.inp";
  print "$cmd\n";
  system "$cmd";
  }
