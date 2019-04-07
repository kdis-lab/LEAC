#!/usr/bin/perl -w
use strict;
use warnings;

my %vardefault =
    (  
       ":SED", 1.0e+16
       , ":SSE", 1.0e+16
       , ":Distortion", 1.0e+16 
       , ":J1",1.0e+16
       , ":Index I", -1000.0 # 0.0
       , ":Index I c-Fuzzy", -1000.0 # 0.0
       , ":CS measure", 1.0e+16
       , ":Dunn's index", 1.0e+16
       , ":SDunn's index", 1.0e+16
       , ":Silhouette", -1000.0 # 0.0
       , ":Simplified Silhouette", -1000.0 # 0.0
       , ":DB-index", 1.0e+16
       , ":Variance Ratio Criterion",  -1000.0 # 0.0
       , ":WB-index",  1.0e+16
       , ":SSB", 1.0e+16
       , ":Score Function", -1000.0 # 0.0
       , ":Rand Index", -1000.0 # 0.0
       , ":Purity", -1000.0 # 0.0
       , ":F-measure", -1000.0 # 0.0
       , ":Jaccard Index", -1000.0 # 0.0
       , ":precision", -1000.0 # 0.0
       , ":recall", -1000.0 # 0.0
       , ":Misclassified", 1.0e+16
       , ":Jm", 1.0e+16
       , ":Xie-Beni index", 1.0e+16
       , ":Xie-Beni index crisp", 1.0e+16
       , ":Entropy", 1000.0 # 10.0
       , "_iterations need", 1.0e+16
       , "_execution time seconds", 1.0e+16
       , "_time seconds need best", 1.0e+16
       , "percent_run_geqk2_tst", 1.0e+16
       , ":DB-index_kge2", 1.0e+16       
    );

my %varcolumn =
    (  
       ":SED", 0
       , ":SSE", 0 
       , ":Distortion", 0
       , ":J1", 0
       , ":Index I", 0
       , ":Index I c-Fuzzy", 0
       , ":CS measure", 0
       , ":Dunn's index", 0
       , ":SDunn's index", 0
       , ":Silhouette", 0
       , ":Simplified Silhouette", 0
       , ":DB-index", 0
       , ":Variance Ratio Criterion", 0
       , ":WB-index",  0
       , ":SSB", 0
       , ":Score Function", 0
       , ":Rand Index", 0
       , ":Purity", 0
       , ":F-measure", 0
       , ":Jaccard Index", 0
       , ":precision", 0
       , ":recall", 0
       , ":Misclassified", 0
       , ":Jm", 0
       , ":Xie-Beni index", 0
       , ":Xie-Beni index crisp", 0
       , ":Entropy", 0
       , "_iterations need", 0
       , "_execution time seconds", 0
       , "_time seconds need best", 0
       , "percent_run_geqk2_tst", 0
       , ":DB-index_kge2", 0
    );


open (FILE,$ARGV[0]) or die "cannot open $ARGV[0] for reading: $!";

my $header = <FILE>;
chomp($header);
my @var = split(/,/,$header);


my $idataset = 0;
my $ialgorithmo = 0;

my $i = 0; 
while($i <= $#var)
{
    if ( ("_datasetkey" eq $var[$i]) or ("_dataset" eq $var[$i])  ) {
	$idataset = $i;
    }
    if ( "_algorithmo" eq $var[$i] ) {
	$ialgorithmo = $i;
    }
    if ( exists $varcolumn{$var[$i]} ) {
	$varcolumn{$var[$i]}  = $i;
    }
    
    $i++;
}
#print "Metrica: $ARGV[0]\n";
# print "$idataset $var[$idataset]\n";
# print "$ialgorithmo $var[$ialgorithmo]\n";
# foreach my $metric (sort keys %varcolumn) {
#     print " $metric $varcolumn{$metric}\n";
# }

my %dataset = ();

my $line;
while ($line = <FILE>) {
	chomp($line);
	my @values = split(/,/,$line);
	$dataset{$values[$idataset]} = ""; 
}
close (FILE);

# foreach my $db (sort keys %dataset) {
#     print " $db $dataset{$db}\n";
# }

foreach my $metric (sort keys %vardefault ) {
    my %totaldatasetvalues = %dataset;
    my $header = "DATASET";
    my $namealgorithmo;
    foreach my $arg (@ARGV) {
	my %datasetvalues = %dataset;
	foreach my $dbvalue (sort keys %datasetvalues ) {
	    $datasetvalues{$dbvalue} = $vardefault{$metric};
	    #print " $db $dataset{$db}\n";
	}	
	open (FILE,$arg) or die "cannot open $arg for reading: $!";
	$line = <FILE>; #LEER HEDER
	while ($line = <FILE>) {
	    chomp($line);
	    my @values = split(/,/,$line);
	    $datasetvalues{$values[$idataset]} = $values[$varcolumn{$metric}];
	    $namealgorithmo = $values[$ialgorithmo];
	}
	$header .= "," . $namealgorithmo;
	close (FILE);
	my $minmax = -1 * $vardefault{$metric} / abs($vardefault{$metric});
	foreach my $dbvalue (sort keys %totaldatasetvalues) {
	    #my $valueout = sprintf("%.3f", $datasetvalues{$dbvalue});
	    if ( $totaldatasetvalues{$dbvalue} eq "" ) {
		$totaldatasetvalues{$dbvalue} = $minmax * $datasetvalues{$dbvalue};
	    }
	    else {
		$totaldatasetvalues{$dbvalue} .= "," . ($minmax * $datasetvalues{$dbvalue});
	    }
	}
    }
    print "------------------------------- $metric $varcolumn{$metric}\n";
    my $namefile = $metric;
    my $letter = substr $metric, 0, 1;
    if ($letter eq "_" or $letter eq ":") {
	$namefile =~ s/^.//;
    }
    $namefile =~ tr/ '-/___/;
    $namefile .= "_table.csv"; 
    $namefile = lc $namefile;
    print "$namefile\n";
    open(my $FILEOUT, '>', $namefile) or die "Could not open file '$namefile' $!";
    print $FILEOUT  "$header\n";
    foreach my $dbvalue (sort keys %totaldatasetvalues ) {
	print $FILEOUT "$dbvalue,$totaldatasetvalues{$dbvalue}\n";
    }
    close($FILEOUT);
}
# foreach $dataset (sort keys %datasetavg) {
#     print "$dataset: $datasetavg{$dataset}\n";
# }
