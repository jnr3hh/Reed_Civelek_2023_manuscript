#! /usr/local/bin/perl

# count links in a series of graphs

#the node file
$in = $ARGV[0];

#the pattern of the series of graphs
$pat = $ARGV[1];
$start = $ARGV[2];
$end =  $ARGV[3];
$output1 = $ARGV[4];
$output2 = $ARGV[5];
$cutoff = $ARGV[6];

srand(1000);
#get list of nodes
$seedF = 1;
$n=0;
open(IN, $in) || die;
while(<IN>) {
    @tmp = split(/[ \t]+/);
    $node{$tmp[0]}=$n;
    $Node[$n] = $tmp[0];
    $n++;
}
close(IN);

#initial count
for($i=0; $i<$n; $i++) {
    for($j=0; $j<$n; $j++) {
	$count[$i][$j] =0.0;
    }
}
print "initiation is done....\n";

#check the series of graphs
$total = 0;
for($i= $start; $i<=$end; $i++) {
    $file = "$pat.$i";
    #print "$file\n";
    if( -s "$file" >0) {
	$total++;
	open(FILE, $file) || die "$file is not found.";
	while(<FILE>) {
	    if(/->/) {
		chop($_);
		#if($seed{$_}>0 || $seedF==0) {
		#    $count{$_}++;
		#}
		@tmp=split('->');
		if(length($node{$tmp[0]}) ==0 || length($node{$tmp[1]}) ==0 ) {
		    print "$_\n";
		    exit(0);
		}
		$count[$node{$tmp[0]}][$node{$tmp[1]}] ++;
	    }
	}
	close(FILE);
    }
}
#clear 
%node =();
#output
print "outputing links...\n";
$hcutoff = $cutoff/2;
open(OUT1, ">$output1") || die;
print OUT1 "digraph G {\n";
for($i=0; $i<$n; $i++) {
    #print "$i\n";
    for($j=0; $j<$n; $j++) {
	$r = $count[$i][$j] /$total;
	$R = $count[$j][$i] /$total;
	if($cutoff>0) {
	    if($r>=$hcutoff && ($R+$r)>=$cutoff && $r>=$R) {
		if($r==$R) {
		    $r += (rand(100)/10000);
		}
		printf(OUT1  "%s->%s [label=%f];\n", $Node[$i],$Node[$j],$r);
	    }
	}
    }
}
print OUT1 "}\n";
close(OUT1);

print "outputing link matrix...\n";
open(OUT2, ">$output2") || die;
for($i=0; $i<$n; $i++) {
    #print "$i\n";
    for($j=0; $j<$n; $j++) {
	$r = $count[$i][$j] /$total;
	$R = $count[$j][$i] /$total;
	if($cutoff>0) {
	    if($r>=$hcutoff && ($R+$r)>=$cutoff && $r>=$R) {
		if($r==$R) {
		    $r += (rand(100)/10000);
		}
		printf( OUT2 "%f\t", $r);
	    }
	    else {
		print OUT2 "0\t";
	    }
	}
	else {
	    printf( OUT2 "%f\t", $r);
	}
    }
    print OUT2 "\n";
}
close(OUT2);

print STDERR "Total $total\n";

