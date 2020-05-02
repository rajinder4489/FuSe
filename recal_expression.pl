#Recalculates the expresion using the SFPG data and the RNA-Seq expression data of transcripts

use strict;
use warnings;
use Data::Dumper;
use Storable;
use List::MoreUtils 'pairwise';
use Getopt::Long;
use Cwd;

################Flags START######################

my ($ori_file, $recal_type, $sfpg, $out_path, $recal_file, $help);

$recal_type = 2;
$recal_file = "recal_exp.txt";
$out_path = getcwd;

GetOptions(
	'input|in=s' => \$ori_file,
	'type|t=i' => \$recal_type,
	'sfpg|g=s' => \$sfpg,
	'out_path|o=s' => \$out_path,
	'recal|re=s' => \$recal_file,
	'help|h' => \$help,
) or die usage();

if($help){usage();}
else
{
	if($recal_type > 2 || $recal_type <= 0)
	{print STDERR "Error Argument '--type' or '-t' should be 1 or 2\n";
		usage();}
	elsif (not defined $ori_file)
	{print STDERR "Error Argument '--input' or '-in' is mandatory\n";
		usage();}
	elsif (not defined $sfpg)
	{print STDERR "Error Argument '--sfpg' or '-g' is mandatory\n";
		usage();}
}

sub usage {
    $0 =~ s/.+\///g;
	print "\nUsage: $0 --input <filename> --sfpg <filename> [--type <1|2>] [--out_path <out path>] [--recal <filename>]\n";
    print "Usage: $0 -in <filename> -g <filename> [-t <1|2>] [-o <out path>] [-re <filename>]\n";
    print "       $0 --help\n";
    print "       $0 -h\n";
	print "
	input|in        Path and filename to the original FPKM file
	
	type|t          Type of recalcualtion
	                Equal distributioin (ED) or Group size distribution (GD)
	                Options 1 for ED or 2 for GD
	                Default 2
	
	sfpg|g          Path and filename to the SFPG data object created using make_sfpgs.pl
	                Default location and name: 'current_folder/int_bla_FuSe_do.data'
	
	out_path|o      Out path for the recalculated FPKM files to be formed
	                Default path is the current folder
	
	recal|re        Filename for the SFPG FPKM out file created using this prog
	                Default location and name: 'current_folder/SFPG_recal_FPKM.txt'
	                An additional file is created with the prefix 'extra' containing the before and after expression of the transcripts that form SFPGs
	
	help|h          To invoke this help\n";
	exit
}

###########################Flags END##################

my $extra_recal_file = "extra_".$recal_file;

print "\nInput\n\tFPKM-> $ori_file\n\tSFPG-> $sfpg\n\tType-> $recal_type\nOutput\n\tRecal SFPGs FPKM-> $recal_file\n\tExtra Recal SFPGs FPKM-> $extra_recal_file\nAt\n\t$out_path\n\n";
##############

print "Reading the expression file\n";
open(my $reads, $ori_file) or die "Could not open file $ori_file $!";

my ($header, %reads);
my $i =0;

while(<$reads>)
{
	$i++;
	chomp;
	if($i == 1)
	{
		$header = $_;
	}
	else
	{
		my @a_r = split /\t/, $_;
		s{\s|\"|\t|\r}{}g foreach @a_r;
		my $id = shift @a_r;
		$reads{$id} = \@a_r;
	}
}


print "\tDone\nLoading SFPG data\n";
my %groups = %{ retrieve("$sfpg") };   # direct to hash
print "\tDone\nCalculating SFPG expression\n";

open(my $fh1, '>', $out_path."/".$recal_file) or die "Could not open file '$recal_file' $!";

open(my $fh2, '>', $out_path."/".$extra_recal_file) or die "Could not open file '$extra_recal_file' $!";

print $fh1 "$header\n";
print $fh2 "type\t$header\n";	#to print the exp of the transcripts before and after the SFPG calculation

my %reads_new;

foreach my $k (sort {lc $a cmp lc $b} keys %groups)
{
	if(exists $reads{$k})
	{
		my @members = split /,/, @{$groups{$k}}[2];
		my $no_mem_main_group = @{$groups{$k}}[1];
		
		my @len = @{$reads{$k}};
		my $no_samples = scalar @len;
		my @sum = ("0") x $no_samples;
		
		foreach my $m (@members)
		{
			$m =~ s/\s//g;
			my $no_groups = @{$groups{$m}}[1];	#number of members in the group of a member; required to calculate the division of reads between the groups;
					
			if(exists $reads{$m})
			{
				my @arr = @{$reads{$m}};
				
				if($recal_type == 1)		#equal distribution
				{
					foreach my $x (@arr) { $x = ($x/$no_groups); }
				}
				elsif($recal_type == 2)		#Number of members distribution
				{
					
					my @mem = split /,/, @{$groups{$m}}[2];	#members in group of a member
					s{\s|\t|\r}{}g foreach @mem;
					
					my $all_num_mem = 0;			#it provides the total number of members for all the groups in which the current transcript id is present
					
					foreach my $mm (@mem)
					{
						$all_num_mem += @{$groups{$mm}}[1];
					}
					
					foreach my $x (@arr) { $x = ($x/$all_num_mem)*$no_mem_main_group; }
						
				}
				else
				{
					print "Wrong recal type\n$recal_type\n";
				}
					
				@sum = pairwise { $a + $b } @sum, @arr;
			}
		}
		$reads_new{$k} = [@sum];
#		die;
	}
}

foreach my $keys (keys %reads)
{
	if(exists $reads_new{$keys})
	{
		print $fh1 "\n$keys";
		
		print $fh2 "\nOriginal\t$keys";
		foreach(@{$reads{$keys}}){print $fh2 "\t$_";}
		
		print $fh2 "\nSFPG\t$keys";
		
		foreach(@{$reads_new{$keys}})
		{
			print $fh1 "\t$_";
			print $fh2 "\t$_";
		}
	}
	else		#the transcripts which do not have SFPGs
	{
		$reads_new{$keys} = @{$reads{$keys}};
		print $fh1 "\n$keys";
		
		foreach(@{$reads{$keys}})
		{
			print $fh1 "\t$_";
		}
	}
}
close $fh1;
close $fh2;

print "\tDone\nFinished\n";
