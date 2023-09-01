#! /usr/bin/perl

$workDir = $ENV{PWD};

$dim = $ARGV[0];
$nNodes = $ARGV[1];
$init = $ARGV[2];
$data = $ARGV[3];
$label = $ARGV[4];

if (length($label)==0) {
    $label ="Label";
}

$prog = "/home/jnr3hh/BN2Distribute/BN/testBN -f 0 ";
$N=1000;

#generate consensus
$cmd="perl countDirectLinksMatrix.pl $data result.out 0 $N junk.0 dag.0 0";
print "$cmd\n";
system($cmd);
unlink("junk.0");
$cutoff=3;
    if (!-e "result.links.$cutoff") {
	$cmd="perl countDirectLinksMatrix.pl $data result.out 0 $N result.links.$cutoff result.linksMatrix.$cutoff 0.$cutoff";
	print "$cmd\n";
	system($cmd);
   
	$cmd="$prog -b $init -d $data -t 0  -D $dim -o result.links$cutoff -c result.links.$cutoff >junk$cutoff";
	print "$cmd\n";
	system($cmd);
    }

