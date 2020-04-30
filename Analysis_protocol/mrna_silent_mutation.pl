#!/usr/bin/perl

use Getopt::Long;
use Data::Dumper;
use POSIX;

my ($help, $mRNA_file, $codon_trans_time_file, $type);
GetOptions(
 'help|h!' => \$help,
 'sequence|s=s' => \$mRNA_file,
 'transtime|t=s' => \$codon_trans_time_file,
 'type|y=s' => \$type,
);

my $usage = "
  Usage: perl mrna_silent_mutation.pl 
              --sequence | -s <mRNA_seq.txt> for mRNA sequence
              --transtime | -t <Trans_time.dat> for individual codon translation time
              --type | -y <mutation type> can be 'slow', 'fast' and 'random'
              [--help | -h]\n\n";

if($help)
{
  die($usage);
}
elsif(!(defined($mRNA_file) && defined($codon_trans_time_file) & defined($type)))
{
  die($usage);
}

my @str = split(/\//, $mRNA_file);
@str = split(/\./, $str[$#str]);
my $file_name = $str[0];
if($type eq "slow")
{
  $type = 1;
  $file_name .= "_slow.txt";
}
elsif($type eq "fast")
{
  $type = 2;
  $file_name .= "_fast.txt";
}
elsif($type eq "random")
{
  $type = 3;
  $file_name .= "_random.txt";
}
else
{
  print("Error: Wrong mutation type.\n");
  die($usage);
}

%codon_table = ("A" => ["GCU", "GCC", "GCA", "GCG"],
				"C" => ["UGU", "UGC"],
				"D" => ["GAU", "GAC"],
				"E" => ["GAA", "GAG"],
				"F" => ["UUU", "UUC"],
				"G" => ["GGU", "GGC", "GGA", "GGG"],
				"H" => ["CAU", "CAC"],
				"I" => ["AUU", "AUC", "AUA"],
				"K" => ["AAA", "AAG"],
				"L" => ["UUA", "UUG", "CUU", "CUC", "CUA", "CUG"],
				"M" => ["AUG"],
				"N" => ["AAU", "AAC"],
				"P" => ["CCU", "CCC", "CCA", "CCG"],
				"Q" => ["CAA", "CAG"],
				"R" => ["AGA", "AGG", "CGU", "CGC", "CGA", "CGG"],
				"S" => ["UCU", "UCC", "UCA", "UCG", "AGU", "AGC"],
				"T" => ["ACU", "ACC", "ACA", "ACG"],
				"V" => ["GUU", "GUC", "GUA", "GUG"],
				"W" => ["UGG"],
				"Y" => ["UAU", "UAC"],
				"Z" => ["UAA", "UAG", "UGA"]);

open(DAT, "<$mRNA_file")||die("Error: cannot find $mRNA_file\n\n");
my $mRNA_seq = "";
while(my $line = <DAT>)
{
  chomp($line);
  $line =~ s/^\s+|\s+$//g;
  $mRNA_seq .= $line;
}
close(DAT);

if(length($mRNA_seq) % 3 ne 0)
{
  die("Error: mRNA sequence is not the multiply of 3\n\n");
}

my $codon_num = length($mRNA_seq) / 3;
print "mRNA sequence has $codon_num codons\n";

my @codon_list = ();
for(my $i = 0; $i < $codon_num; $i++)
{
  push@codon_list, substr($mRNA_seq, 3*$i, 3);
}
my @aa_list = ();
foreach my $c (@codon_list)
{
  push@aa_list, codon2AA($c);
}

open(DAT, "<$codon_trans_time_file")||die("Error: cannot find $codon_trans_time_file\n\n");
my %ctt_table = ();
while(my $line = <DAT>)
{
  chomp($line);
  my @str = split(/\s+/, $line);
  my $codon = $str[0];
  my $ctt = $str[1];
  $codon =~ s/^\s+|\s+$//g;
  $ctt =~ s/^\s+|\s+$//g;
  $ctt_table{$codon} = $ctt;
}
close(DAT);

@new_codon = ();
foreach my $a (@aa_list)
{
  my @c_list = @{$codon_table{$a}};
  if($type eq 1) # slow
  {
  	my $t = 0;
  	my $cc = "";
  	foreach my $c (@c_list)
  	{
  	  if($ctt_table{$c} > $t)
  	  {
  	  	$t = $ctt_table{$c};
  	  	$cc = $c;
  	  }
  	}
  	push@new_codon, $cc;
  }
  elsif($type eq 2) # fast
  {
  	my $t = 99999999;
  	my $cc = "";
  	foreach my $c (@c_list)
  	{
  	  if($ctt_table{$c} < $t)
  	  {
  	  	$t = $ctt_table{$c};
  	  	$cc = $c;
  	  }
  	}
  	push@new_codon, $cc;
  }
  elsif($type eq 3) # random
  {
  	my $idx = int(rand(@c_list));
  	push@new_codon, $c_list[$idx];
  }
}

open(DAT, ">$file_name")||die("Error: cannot create $file_name");
for(my $i = 0; $i < POSIX::ceil(@new_codon/20); $i++)
{
  my $b = (($i+1)*20, $#new_codon+1)[($i+1)*20 > @new_codon];
  for(my $j = $i*20; $j < $b; $j++)
  {
  	print DAT $new_codon[$j];
  }
  print DAT "\n";
}

########################################################################
sub codon2AA
{
  my $codon = $_[0];
  for my $a (keys(%codon_table))
  {
  	if($codon ~~ @{$codon_table{$a}})
  	{
  	  return $a;
  	}
  }
}