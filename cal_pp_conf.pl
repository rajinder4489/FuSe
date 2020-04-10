#Rajinder Gupta
#This program finds the confidence of protein pairs
#It only keeps the cases where both the proteins in the pair are protein coding

use strict;
use warnings;
use Data::Dumper;
#use Array::Utils qw(:all);
use Storable;
use List::Util qw(sum0);
use Getopt::Long;
use Cwd;

################Passing flags START######################

my ($rel_file, $bi_do, $ss, $out_path, $pp_do, $help);

$rel_file = "data/Biomart_rel.txt";
$bi_do = "data/bi_do.data";
$ss = "data/scoring_scheme.txt";
$pp_do = "pp_FuSe.data";
$out_path = getcwd;

GetOptions(
	'rel|r=s' => \$rel_file,
	'bi_do|bi=s' => \$bi_do,
	'ss|s=s' => \$ss,
	'pp_do|p=s' => \$pp_do,
	'out_path|o=s' => \$out_path,
	'help|h' => \$help,
) or die usage();

if($help){usage();}

sub usage {
    $0 =~ s/.+\///g;
	print "\nUsage: $0 [--rel <filename>] [--bi_do <filename>] [--ss <filename>] [--out_path <out path>] [--pp_do <filename>]\n";
    print "Usage: $0 [-r <filename>] [-bi <filename>] [-s <filename>] [-o <out path>] [-p <filename>]\n";
    print "       $0 --help\n";
    print "       $0 -h\n";
	print "
	rel|r           Path and filename to the Biomart relation file
	                Default location and name: 'data/Biomart_rel.txt'
	
	bi_do|bi        Path and filename to the BLAST and Interpro data object created using preprocessing_data.pl
	                Default location and name: 'data/bi_do.data'
	
	ss|s            Path and filename to the scoring scheme file
	                Default location and name: 'data/scoring_scheme.txt'
	
	pp_do|o         Filename for the data object to be created using this program
	                Default name 'pp_FuSe.data'
	
	out_path|o      Out path for the pp_FuSe.data
	                Default path is the current folder
	
	help|h          To invoke this help\n";
	exit
}

###########################Passing flags END##################
print "\nInput\n\tBiomart relation file-> $rel_file\n\tBLAST Interpro do-> $bi_do\n\tScoring scheme-> $ss\nOutput\n\tData object-> $pp_do\nAt\n\t$out_path\n\n";

###########

open(my $rel, $rel_file) or die "Could not open file '$rel_file' $!";

my (%pg, %pt, %pty, %c, %pg_combo, %ppt, %blast, %no_match, %match);	#pg- protein-gene, pt- protein-transcript, pty - protein-type, c- protein-confidence, pg_combo - protein pair - gene pair combo, ppt - protein pair type (isozyme/isoform), no_match- interpro db which show no match for protein seq, match- interpro db which show match for protein seq

while(<$rel>)
{
	chomp;
	my @arr = split /\t/, $_;
	s{\n|\r|\t|\s}{}g foreach @arr;
	
	$pg{$arr[0]} = $arr[1];
	$pt{$arr[0]} = $arr[2];
	$pty{$arr[0]} = $arr[3];
}

print "Reading BLAST-Interpro data object\n";
my %hash = %{ retrieve($bi_do) };   # direct to hash		#gets the data object saved from the last script

print "\tDone\n\n\n";

print "Reading scoring scheme\n";
my $file_ss = $ss;
open(my $f_ss, $file_ss) or die "Could not open file '$file_ss' $!";

my @ss_blast;
my %ss_hash;
my %tool_pos;

my $pos = 4;		#The first 4 values (0,1,2,3) in the array with pair pair values are BLAST identity, sub cov, que cov, and gap%.

while(<$f_ss>)
{
	chomp;
	
	if($_ =~ m/\>/)
	{
		my @a = split /\t/, $_;
		s{\s+|\r|\n|\>}{}g foreach @a;
		
		if($_ =~ m/BLAST/i)
		{
			@ss_blast = ($a[2], $a[3], $a[4]);
			$ss_hash{$a[0]} = [$a[2], $a[3], $a[4]];
		}
		else
		{
			$ss_hash{$a[0]} = [$a[1], $a[2], $a[3], $a[4], $a[5]];
#			$tool_pos{$a[0]} = $pos;
			$tool_pos{$pos} = $a[0];
			$pos++;
		}
	}
}
print "\tDone\n\nCalculating protein pair (pp) confidence scores\n";

my %pairs;
my ($zz, $yy, $xx, $zz1, $xx1);		#zz - protein pairs both prot coding, yy - used for printing, xx - protein pairs atleast one not prot coding
$zz = $yy = $xx = $zz1 = $xx1 = 0;

foreach my $k (keys %hash)
{
	my @arr = @{$hash{$k}};
	my @arr1 = @arr;

	s{STONM}{0}g foreach @arr;
	s{STNM}{1}g foreach @arr;
	s{STM}{2}g foreach @arr;
	s{NM}{3}g foreach @arr;
	s{NP12|NP1|NP2}{4}g foreach @arr;
	
#	Score1 is knowledege limited (seq and domains); for development only
	my $score1 = ($arr[0]*$ss_blast[0]/100) +	#blast identity
		((100 - $arr[1])*$ss_blast[1]/100) +	#blast coverage for sequence 1 (query seq)
		((100 - $arr[2])*$ss_blast[1]/100) + 	#blast coverage for sequence 2 (subject seq)
		($arr[3]*$ss_blast[2]/10);				#gapopen is added to the formula
	
#	Score2 is knowledge score; normalizes for the cases where no information is present from Interpro tool(s)
	my $score2 = $score1;
	
#	Score3 discovery score; it is alignment dependent only
	my $score3 = $score1 * 100/$ss_blast[0];				#The last part normalizes the score to a scale of 100
	
	my $max = $ss_blast[0];
	for(my $i=4; $i<=$#arr; $i++)							#for score1 and score2
	{
		$score1 += @{$ss_hash{$tool_pos{$i}}}[$arr[$i]];		#$tool_pos{$i} -> gives the name of the tool; this name is then used by $ss_hash to get to the ss for that tool; $arr[$i] gives if it is STONM, STOM .... etc.
		
		if($arr[$i] != 4)
		{
			$score2 += @{$ss_hash{$tool_pos{$i}}}[$arr[$i]];
			$max += @{$ss_hash{$tool_pos{$i}}}[0];
		}
	}

	$score2 = ($score2/$max)*100;
	
	my @p = split /\-/, $k;
	s{\.[0-9]+}{}g foreach @p;
	
	my $t1 = $pt{$p[0]};
	my $t2 = $pt{$p[1]};
	
	$zz++;
	if($zz>=$yy+1000000)
	{
		print "\tCalculated for $zz pp\n";
		$yy = $zz;
	}
	
	if($pty{$p[0]} eq 'protein_coding' && $pty{$p[1]} eq 'protein_coding')	#removes all other type of sequences
	{
		$pairs{$k} = [$t1, $t2, $score1, $score2, $score3, @arr1];	
	}
}
print "\tCalculated for $zz\nDone\n\n";

print "Saving the data object $pp_do\n";
store \%pairs, join '/', $out_path, $pp_do or die;

print "\tDone\nFinished\n";