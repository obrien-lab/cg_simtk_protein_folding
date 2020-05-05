#!/usr/bin/perl

$t=$ARGV[0];
$p=$ARGV[1];
$kbT=1.9872/1000*$t;

$dGns = -$kbT*log($p/(1-$p));
print "$dGns kcal/mol at $t K\n";
