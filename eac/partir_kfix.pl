#!/usr/bin/perl -w
open (FILE,$ARGV[0]) or die "cannot open $ARGV[0] for reading: $!";
while (<FILE>) {
    #chomp
    @var = split(",");
    @dataset = split("/",$var[0]);
    print "_datasetkey,$dataset[6]";
#    print "_datasetkey,$dataset[4]_$dataset[6]";
    for ($i=3; $i < @var; $i++) {
	print ",$var[$i]" 
	    #    $char = substr($var[$i], 0, 1);
	    #    if ($char eq "_" or $char eq ":") {
	    #	$iVar = $var[$i];
	    #	$iValue = $var[$i+1];
	    #	print "$i\t$iVar = $iValue\n";
	    #    }
    }
}
close(FILE);

