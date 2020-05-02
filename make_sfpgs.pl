#makes the SFPG based on csc and protein pair data object 

use strict;
use warnings;
use Storable;
use Data::Dumper;
use List::Util qw(sum);
use Getopt::Long;
use Cwd;

################Flags START######################

my ($pp_do, $conf, $score_type, $out_path, $sfpg_file, $help);

$conf = 95;
$out_path = getcwd;
$score_type = "KS";

GetOptions(
	'pp_do|p=s' => \$pp_do,
	'score_type|t=s' => \$score_type,
	'csc|c=i' => \$conf,
	'sfpg|g=s' => \$sfpg_file,
	'out_path|o=s' => \$out_path,
	'help|h' => \$help,
) or die usage();

my %st = ('KS', 3, 'ks', 3, 'Ks', 3,  'kS', 3, 'DS', 4, 'ds', 4, 'Ds', 4,  'dS', 4, 'KDS', 2);

if($help){usage();}
else
{
	if(not defined $pp_do)
	{print STDERR "Error Argument '--pp_do' or '-p' is mandatory\n";
		usage();}
	elsif($conf<=0 || $conf>100)
	{print STDERR "Error Argument '--conf' or '-si' should be >0 & ≤100\n";
		usage();}
	elsif(!exists $st{$score_type})
	{print STDERR "Error Valid values for argument '--score_type' or '-t' are 'KS' or 'DS'\nBoth upper and lower case are allowed\n";
		usage();}
}

sub usage {
    $0 =~ s/.+\///g;
	print "\nUsage: $0 --pp_do <filename> [--score_type <KS|DS>] [--conf <integer>] [--out_path <out path>] [--sfpg <filename>]\n";
    print "Usage: $0 -p <filename> [-t <KS|DS>] [-c <integer>] [-o <out path>] [-g <filename>]\n";
    print "       $0 --help\n";
    print "       $0 -h\n";
	print "
	pp_do|o         Path and filename to the data object created using cal_pp_conf.pl
	                	
	score_type|t    The type of score to use for calculation of FPKM for SFPGs
	                Options KS | DS
	                Both upper and lower case allowed (KS | kS | Ks | ks | DS | dS | Ds | ds)
					Default 'KS'
	
	conf|c          Confidence score cutoff, to take only the protein pairs above it for forming SFPGs
	                Positive Integer; range >0 to ≤100
	                Default ≥95
	
	sfpg|g          Filename for the data object to be created using this program
	                Default value 'sfpg_<KS|DS>_<conf>.data'
					An additional file with the same name is created as .txt
	
	out_path|o      Out path for the SFPG data object to be formed
	                Default path is the current folder
	
	help|h          To invoke this help\n";
	exit
}

###########################Flags END##################
my ($out_file1, $out_file2);

if(not defined $sfpg_file)
{
	$out_file1 = "sfpg_".$score_type."_".$conf.".txt";
	$out_file2 = "sfpg_".$score_type."_".$conf.".data";
}
else
{
	if($sfpg_file=~ m/\.data/){}
	else{$sfpg_file = $sfpg_file.".data"}
	
	$out_file1 = $sfpg_file;
	$out_file2 = $sfpg_file;
	$out_file1 =~ s/.data/.txt/;
}

print "\nInput\n\tProtein pairs-> $pp_do\n\tSimilarity score cutoff-> $conf\nOutput\n\tSFPGs-> $out_file2\n\t        $out_file1\nAt\n\t$out_path\n\n";
##############

print "Loading the Protein pair data\n";
my %pairs = %{ retrieve("$pp_do") };   # direct to hash
print "\tDone\n";

#my $out_file1 = $out_path.$out_file1;
open(my $out_fh1, '>', $out_path."/".$out_file1) or die "Could not open file '$out_file1' $!";

my (%groups, %avscore);		#average score for the group

my $ti = $st{$score_type};

print "Making SFPGs\n";

foreach my $k (keys %pairs)
{
	my @a = @{$pairs{$k}};
	s{\n|\r|\s}{}g foreach @a;

	if($a[$ti]>= $conf)
	{
		if(exists $groups{$a[0]})
		{
			push @{$groups{$a[0]}}, $a[1];
			push @{$avscore{$a[0]}}, $a[$ti];
		}
		else
		{
			$groups{$a[0]} = [$a[1]];
			$avscore{$a[0]} = [$a[$ti]];
		}
		
		if(exists $groups{$a[1]})
		{
			push @{$groups{$a[1]}}, $a[0];
			push @{$avscore{$a[1]}}, $a[$ti];
		}
		else
		{
			$groups{$a[1]} = [$a[0]];
			$avscore{$a[1]} = [$a[$ti]];
		}
	}
}
print "\tDone\nCalculating SFPG confidence\n";

my %g_formed;
foreach my $k (keys %groups)
{
	my $c = scalar @{$avscore{$k}}+1;
	my $av = (sum @{$avscore{$k}})/(scalar @{$avscore{$k}});
	my $str = $k.', '.join ', ', @{$groups{$k}};
	$g_formed{$k} = [$k, $c, $str, $av];
	print $out_fh1 "$score_type\t$k\t$c\t$str\t$av\n";
}

print "\tDone\n\nSaving the data object\n";
store \%g_formed, $out_path."/".$out_file2 or die;

print "\tDone\nFinished\n";
