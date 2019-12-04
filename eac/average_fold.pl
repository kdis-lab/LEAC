#!/usr/bin/perl -w

#use warnings;
#use strict;
#use Scalar::Utils
#    qw(looks_like_numbers);


# 1	_inout = out
# 7	_times run = 1
# 13	_random seed = 3874402244 1419207497 4133933127 3971050141 822347760 4199991182 1147493713 3616605934
# 15	_maximum execution time = 36000
# 17	_print centroids format = false
# 19	_print table format = false
# 23	, "_number maximum generations = 100, 1
# 25	_probability crossover = 0.8
# 27	_probability mutation = 0.001
# 33	_format file = keel
# 41	_select-attributes = 1-8, 1
# 43	_class-column = 10, 1
# 45	:dataset test = /home/user/dataset/abalone6cf_ok/abalone6cf-10-10tst.dat
# 49	_number run = 1
# 53	_end condition = true
# 55	_end condition test = true
# 57	_run time message = ok
# 59	_file name plot funtion hist = none


$varsrt = "";
$countClass = 0;

# 1 SE PUEDE CALCULAR INDEPENDIENTE DE UNA INDETERMINACION
# 2 EN CASO DE UNA INDETERMINACION SE LE ASIGNAN EL VALOR POR DEFAULT
     #   "_algorithmo", 1
     # , "_author",  1
     # , "_datasetkey", 1
     # , "_d", 1
     # , "_norm", 1
     # , "_size population", 1 
     # , "_number maximum generations", 1 
     # , "_k-minimum", 1
     # , "_k-maximum", 1
     # , "_metric used", 1

%varprocUni = 
    ("_algorithmo",  1
     , "_author", 1
     , "_runnig date", 1
     , "_norm", 1
     , "_size population", 1
     , "_k-minimum", 1
     , "_k-maximum", 1
     , "_datasetkey", 1
     , "_d", 1
     , "_metric used", 1
#     , "_outKClass",1
    );

%varproc = 
    (  
       # "_execution time seconds", 1 
       #  , "_iterations need", 1
       #  , "_time seconds need best", 1
       #  , "_outK", 1
       #  , ":outK", 1
       #  , "_solution overridden in the run of the algorithm", 1 
       #  , "_objetivefuncrun", 2
       #  , "_fitness", 2
       #  , "_number total generations", 1 
       #  , "_n", 1

       "_execution time seconds", 1
       , "_iterations need", 1
       , "_time seconds need best", 1
       , "_outK", 1
       , ":outK", 1
       , "_solution overridden in the run of the algorithm", 1
       , "_objetivefuncrun", 2
       , "_fitness", 2
       , "_number total generations", 1
       , "_n", 1
       
       # , "_SED", 2
       # , "_SSE", 2
       # , "_Distortion", 2
       # , "_J1", 2
       # , "_Index I", 2
       # , "_CS measure", 2
       # , "_Dunn's index", 2
       # , "_SDunn's index", 2
       # , "_Silhouette", 2
       # , "_Simplified Silhouette", 2 
       # , "_DB-index", 2
       # , "_Variance Ratio Criterion", 2 
       # , "_WB-index", 2
       # , "_SSB", 2
       # , "_Score Function", 2
       # , "_Rand Index", 2
       # , "_Purity", 2
       # , "_F-measure", 2
       # , "_Jaccard Index", 2
       # , "_precision", 2
       # , "_recall", 2
       # , "_Misclassified", 2 #, "_Error Bezdek", 1
       # , "_pairs a", 2
       # , "_pairs b", 2
       # , "_pairs c", 2
       # , "_pairs d", 2
       # , "_Jm", 2
       # , ":Rand Index", 2
       
       , "_SED", 2
       , "_SSE", 2
       , "_Distortion", 2
       , "_J1", 2
       , "_Index I", 2
       , "_Index I c-Fuzzy", 2

       , "_CS measure", 2
       , "_Dunn's index", 2
       , "_SDunn's index", 2
       , "_Silhouette", 2
       , "_Simplified Silhouette", 2
       , "_DB-index", 2
       , "_Variance Ratio Criterion", 2

       , "_WB-index", 2
       , "_SSB", 2
       , "_Score Function", 2 
       , "_Rand Index", 2
       , "_Purity", 2
       , "_F-measure", 2,
       , "_Jaccard Index", 2

       , "_precision", 2
       , "_recall", 2
       , "_Misclassified", 2
       , "_pairs a", 2
       , "_pairs b", 2
       , "_pairs c", 2
       , "_pairs d", 2
       , "_Jm", 2
       , "_Xie-Beni index", 2
       , "_Xie-Beni index crisp", 2
       , "_Entropy", 2
       , "_Partition Coefficient", 2
       , "_Intra- and inter-clust", 2
      
       , ":SED", 2
       , ":SSE", 2
	, ":Distortion", 2
	, ":J1", 2
	, ":Index I", 2
        , ":Index I c-Fuzzy", 2
	, ":CS measure", 2
	, ":Dunn's index", 2
	, ":SDunn's index", 2
	, ":Silhouette", 2
	, ":Simplified Silhouette", 2
	, ":DB-index", 2
	, ":Variance Ratio Criterion", 2
	, ":WB-index", 2
	, ":SSB", 2
	, ":Score Function", 2
	, ":Rand Index", 2
	, ":Purity", 2
	, ":F-measure", 2
	, ":Jaccard Index", 2
	, ":precision", 2
	, ":recall", 2
	, ":Misclassified", 2
	, ":pairs a", 2
	, ":pairs b", 2
	, ":pairs c", 2
	, ":pairs d", 2
	, ":Jm", 2
	, ":Xie-Beni index", 2
        , ":Xie-Beni index crisp", 2
	, ":Entropy", 2
	, ":Partition Coefficient", 2
	, ":Intra- and inter-clust", 2
	, ":n", 1       
       
       # , "_Avg Entropy", 2 
       # , "_Xie-Beni index", 2
       # , "_objetivefuncrun", 2
       # , ":SSE", 2
       # , ":TWCV", 2
       # , ":Index I", 2
       # , ":Distortion distance", 2 
       # , ":J1", 2
       # , ":JM", 2
       # , ":Index I", 2
       # , ":Jaccard Index", 2
       # , ":DB-index", 2
       # , ":Simplified Silhouette", 2 
       # , ":Silhouette", 2
       # , ":Avg Purity", 2
       # , ":Avg Entropy", 2
       # , ":Variance Ratio Criterion", 2 
       # , ":CS measure", 2
       # , ":Dunn's index", 2
       # , ":SDunn's index", 2
       # , ":Score Function", 2
       # , ":Xie-Beni index", 2
       # , ":Error Bezdek", 1
       # , ":precision", 2
       # , ":recall", 2
       # , ":F-measure", 2
       # , ":n", 1
    );

%varprocmult =
    (
     "_sensitivity", 2 
     , "_specificity", 2 
     , ":sensitivity", 2 
     , ":specificity" , 2 
     # 141	_sensitivity = E;38.40579710144927;B;35.80092011710581;C;29.00274473924977;F;0;D;0;A;0, 2
     # 143	_specificity = E;14.44141689373297;B;87.88501026694045;C;39.18417799752781;F;0;D;0;A;0, 2
     # 207	:sensitivity = E;34.23423423423424;B;39.51219512195122;C;39.68253968253968;F;0;D;39.47368421052632;A;0, 2
     # 209	:specificity = E;46.91358024691358;B;75;C;27.77777777777778;F;0;D;23.80952380952381;A;0, 2
    );

%vardefault =
    (  
       "_objetivefuncrun", 1.0e+16
       , "_fitness", -1
       , "_SED", 1.0e+16
       , "_SSE", 1.0e+16
       , "_Distortion", 1.0e+16
       , "_J1", 1.0e+16
       , "_Index I", -1000.0 # 0.0
       , "_Index I c-Fuzzy", -1000.0
       # , "_time seconds need best", 1.0e+16

       , "_CS measure", 1.0e+16
       , "_Dunn's index", -1000.0 # 0.0
       , "_SDunn's index", -1000.0 # 0.0
       , "_Silhouette", -1000.0 # 0.0
       , "_Simplified Silhouette", -1000.0 # 0.0
       , "_DB-index", 1.0e+16
       , "_Variance Ratio Criterion", -1000.0 # 0.0

       , "_WB-index", 1.0e+16
       , "_SSB", 1.0e+16
       , "_Score Function", -1000.0 # 0.0
       , "_Rand Index", -1000.0 # 0.0
       , "_Purity", -1000.0 # 0.0
       , "_F-measure", -1000.0 # 0.0
       , "_Jaccard Index", -1000.0 # 0.0
       , "_precision", -1000.0 # 0.0
       , "_recall", -1000.0 # 0.0
       , "_Misclassified", 1.0e+16
       # , "_Error Bezdek", 1.0e+16
       , "_Jm", 1.0e+16
       , "_Xie-Beni index", 1.0e+16
       , "_Xie-Beni index crisp", 1.0e+16
       , "_Entropy", 1000.0 # 10.0
       # , "_time seconds need best", 1.0e+16

       , ":SED", 1.0e+16
       , ":SSE", 1.0e+16
       , ":Distortion", 1.0e+16 
       , ":J1",1.0e+16
       , ":Index I", -1000.0 # 0.0
       , ":Index I c-Fuzzy", -1000.0
       , ":CS measure", 1.0e+16
       , ":Dunn's index", -1.0e+16
       , ":SDunn's index", -1.0e+16
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
    );


foreach $vark (sort keys %varprocUni) {
    $varsrt .= $vark . "," ; 
}
foreach $vark (sort keys %varproc) {
    $varsrt .= $vark . ",," ; 
}
#Percentage runs with k greater than two
$varsrt .= "percent_run_geqk2_tra,percent_run_geqk2_tst,";
$varsrt .= ":DB-index_kge2,,";
$varsrt .=  "_NumClass,";
#foreach $varkm (sort keys %hvarmult) {
foreach $vark (sort keys %varprocmult) {
    $varsrt .= $vark . "," ; 
}    
print "$varsrt\n";
foreach my $arg (@ARGV) {
    open (FILE,$arg) or die "cannot open $arg for reading: $!";
    # $isFist = 0;
    %hvar = ();
    %hvarmult = ();
    $totaltrainingok = 0;
    $totaltestok = 0;
    $db_index = "";
    $totalrun = 0;
    while (<FILE>) {
	chomp;
        # ++$isFist;
	@var = split(/,/);
	%hvarline = ();
	for ($i=0; $i < @var; $i++) {
	    $char = substr($var[$i], 0, 1);
	    if ($char eq "_" or $char eq ":") {
		$hvarline{$var[$i]} = $var[$i+1];
		# print " '$var[$i]'  = '$var[$i+1]'\n";
		$i++;
	    }
	}
	$totalrun++;
	if ( ($hvarline{"_outK"} > 1) and ($hvarline{"_end condition"} eq "true") ) {
	    $totaltrainingok++;
	}
	else {                                               
	    $hvarline{"_time seconds need best"} = $hvarline{"_execution time seconds"}; 
	    $hvarline{"_iterations need"} = $hvarline{"_number maximum generations"};
	    $hvarline{"_Misclassified"} = $hvarline{"_n"};
	    $hvarline{":Misclassified"} = $hvarline{":n"};
	    
	}
	foreach $varline (sort keys %hvarline) {
	    $char = substr($varline, 0, 1);
	    #if (  ($hvarline{"_outK"} > 1) and ($hvarline{"_end condition"} eq "true") ) {
	    if ( $char eq "_" ) {
		#VARIABLES UNICAS
		if ( exists $varprocUni{$varline} ) {
		    if ( !exists $hvar{$varline} ) {
			$hvar{$varline} = $hvarline{$varline};
		    }
		}
		#VARIABLES
		elsif ( (exists $varproc{$varline} and (($hvarline{"_outK"} > 1) and ($hvarline{"_end condition"} eq "true")) ) or (exists $varproc{$varline} and  ($varproc{$varline} eq 1) ) ) {
		    $hvar{$varline} .= $hvarline{$varline} . ",";
		    # print "$varline = $hvar{$varline}\n";
		}
		#VARIABLES MULTIPLES
		elsif ( exists $varprocmult{$varline} and (($hvarline{"_outK"} > 1) and ($hvarline{"_end condition"} eq "true")) ) {
		    @varmulti = split(/;/,$hvarline{$varline});
		    $countClass = 0;
		    for ($ll=0; $ll < @varmulti; $ll+=2) {
			$keymult = $varline ."_" . $varmulti[$ll];
			$hvarmult{$keymult} .= $varmulti[$ll+1] . ",";	    
			#$hvar{$iVar} .= $iValue . ",";
			#print "$ll\t$keymult = $varmulti[$ll+1] = $hvarmult{$keymult}\n";
			$countClass+= 1;
		    }
		    #$varprocUni{"_outKClass"} = $countClass; 
		}
	    } #IF _	
	} #FOREACH
	#} #IF _K > 1

	#if ( ($hvarline{":outK"} > 1) and ($hvarline{"_end condition test"} eq "true") ) {
	if ( ($hvarline{":outK"} > 1) and ($hvarline{"_end condition"} eq "true") ) {
	    $totaltestok++;
	    $db_index .= $hvarline{":DB-index"} . ",";
	}
	foreach $varline (sort keys %hvarline) {
	    $char = substr($varline, 0, 1);
	    if ( $char eq ":" ) {
		#VARIABLES UNICAS
		if ( exists $varprocUni{$varline} ) {
		    if ( !exists $hvar{$varline} ) {
			$hvar{$varline} = $hvarline{$varline};
		    }
		}
		#VARIABLES
		elsif ( (exists $varproc{$varline} and (($hvarline{"_outK"} > 1) and ($hvarline{"_end condition"} eq "true"))) or ( exists $varproc{$varline} and ($varproc{$varline} eq 1)) ) {
		    $hvar{$varline} .= $hvarline{$varline} . ",";
		    # print "$varline = $hvar{$varline}\n";
		}
		#VARIABLES MULTIPLES
		elsif ( exists $varprocmult{$varline} and (($hvarline{"_outK"} > 1) and ($hvarline{"_end condition"} eq "true")) )  {
		    @varmulti = split(/;/,$hvarline{$varline});
		    #$countClass = 0;
		    for ($ll=0; $ll < @varmulti; $ll+=2) {
			$keymult = $varline ."_" . $varmulti[$ll];
			$hvarmult{$keymult} .= $varmulti[$ll+1] . ",";	    
			#$hvar{$iVar} .= $iValue . ",";
			#print "$ll\t$keymult = $varmulti[$ll+1] = $hvarmult{$keymult}\n";
			#$countClass+= 1;
		    }    
		    #$varprocUni{"_outKClass"} = $countClass; 
		}
	    } #IF :	
	} #FOREACH
	#} #IF :K > 1	
    } #WHILE	
    #foreach $vark (sort keys %hvar) {
    #foreach $vark (keys %hvar) {
    $varstat = "";
    #LAS UNICAS
    foreach $vark (sort keys %varprocUni) {
	if (exists $hvar{$vark} ) {
	    $varstat .= $hvar{$vark} . ",";
	}
	else {
	    $varstat .=  ",";
	}
    }
    #ESTADISTICAS
    foreach $vark (sort keys %varproc) {
	$sum = 0;
	$count = 0;
	$devstd = 0;
	if (exists $hvar{$vark} ) {
	    #print "hvar{$vark} = $hvar{$vark}\n";
	    @values = split(/,/, $hvar{$vark});
	    #
	    foreach $iData (@values) {
		#if ( $iData =~ /^[0-9,.E]+$/ ) {
		#if ( $iData =~ /^-?[0-9,.Ee]+$/ ) {
		if ( $iData =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) {
		    $sum += $iData;
		    $count++; 
		}
	    }
	    if ( $count > 0 and $sum =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) {
		$average = $sum / $count;
		
		foreach $iData (@values) {
		    #if ( $iData =~ /^-?[0-9,.Ee]+$/ ) {
		    if ( $iData =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) {
			$difave = ($average - $iData) ;
			$devstd += $difave * $difave; 
		    }
		}
		if ( $count > 1 ) {
		    $devstd = sqrt($devstd / ($count -1));
		}
		else {
		    $devstd = "NaN";
		}
                
	    }
	    else {
		$average = $vardefault{$vark};  # "NaN";
		$devstd  = "NaN";
	    }
	    #print "$vark\t$hvar{$vark} = $average \n";
	    # print "$vark\n";
	    # print "$average\n";
	    #print "$devstd\n";
	    #$varsrt .= $vark . ",";
	    
	    if ($sum =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) {
		#if ( $average > 1.0e+16 ) {
		#    $average = 1.0e+16;
		#}
		if ($average < 1.0)  {
		    $average = sprintf("%.3g", $average);
		}
		else {
		    $average = sprintf("%.3f", $average);
		}
	    }
	    $varstat .= $average . "," . $devstd . ",";
	}
	else {
	    if (exists $vardefault{$vark} ) {
		$varstat .=  $vardefault{$vark} . ",NaN,";
	    }
	    else {
		$varstat .=  "NaN,NaN,";
	    }
	     
	}
    }
    #Percentage runs with k greater than two
    $totaltrainingok = 100 * $totaltrainingok / $totalrun;
    #$totaltrainingok = sprintf("%.3f", $totaltrainingok);
    if ( $totaltrainingok < 1.0 ) {
	$totaltrainingok = sprintf("%.3g", $totaltrainingok);
    }
    else {
	$totaltrainingok = sprintf("%.3f", $totaltrainingok);
    }
    $totaltestok = 100 * $totaltestok / $totalrun;
    if ( $totaltestok < 1.0) {
	$totaltestok = sprintf("%.3g",$totaltestok);
    }
    else {
	$totaltestok = sprintf("%.3f",$totaltestok);
    }
    $varstat .= $totaltrainingok . "," . $totaltestok . ",";

    #ONLY DB-INDEX
    $sum_db = 0;
    $count_db = 0;
    $devstd_db = 0;
    $average_db = 0;
    #print "DB_INDEX = $db_index\n";
    @values_db = split(/,/,$db_index);
    foreach $iData (@values_db) {
	if ( $iData =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) {
	    $sum_db += $iData;
	    $count_db++; 
	}
    }
    if ( $count_db > 0  and  $sum_db =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) {
	$average_db = $sum_db / $count_db;	
	foreach $iData (@values_db) {
	    if ( $iData =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) {
		$difave = ($average_db - $iData) ;
		$devstd_db += $difave * $difave; 
	    }
	}
	if ( $count_db > 1 ) {
	    $devstd_db = sqrt($devstd_db / ($count_db -1));
	}
	else {
	    $devstd_db = "NaN";
	}
	if ( $average_db < 1.0 ) {
	    $average_db = sprintf("%.3g", $average_db);
	}
	else {
	    $average_db = sprintf("%.3f", $average_db);
	}
    }
    else {
	$average_db = $vardefault{":DB-index"};  # "NaN";
	$devstd_db  = "NaN";
    }
    $varstat .= $average_db . "," . $devstd_db . ",";

    #NUM OF CLASS
    $varstat .=  $countClass . ",";;
    # if ( keys(%hvarmult) == 0 ) {
    # 	$varstat .=  "NaN,";
    # }
    # else {
    # 	$varstat .= keys(%hvarmult) / 4  . ",";
    # }
    
    #ESTADISTICAS MULTIPLES
    foreach $varkm (sort keys %hvarmult) {
	$sum = 0;
	$count = 0;
	$devstd = 0;
	#print "STADISTICAS MULTIPLES: $varkm  \n";
	#if (exists $hvar{$vark} ) {
	#print "$varkm = $hvarmult{$varkm}\n";
	@values = split(/,/, $hvarmult{$varkm});
	#
	foreach $iData (@values) {
	    #if ( $iData =~ /^[0-9,.E]+$/ ) {
	    #if ( $iData =~ /^-?[0-9,.Ee]+$/ ) {
	    if ( $iData =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) {
		$sum += $iData;
		$count++; 
	    }
	}
	if ( $count > 0 ) {
	    $average = $sum/ $count;
	    
	    foreach $iData (@values) {
		#if ( $iData =~ /^-?[0-9,.Ee]+$/ ) {
		if ( $iData =~ /^([+-]?)(?=\d|\.\d)\d*(\.\d*)?([Ee]([+-]?\d+))?$/ ) {
		    $difave = ($average - $iData) ;
		    $devstd += $difave * $difave; 
		}
	    }
	    if ($count > 1 ) {
		$devstd = sqrt($devstd / ($count -1));
	    }
	    else {
		$devstd = "NaN";
	    }
	}
	else {
	    $average = "NaN";
	    $devstd  = "NaN";
	}
	#print "$vark\t$hvar{$vark} = $average \n";
	#print "$vark\n";
	#print "$average\n";
	#print "$devstd\n";
	#$varsrt .= $vark . ",";
	$varstat .= $varkm . "," . $average . "," . $devstd . ",";
	#}
	#else {
	#    $varstat .=  ",,";
	#}
    }

    #print "$varsrt\n";
    print "$varstat\n";

    
    #print "N = $n\n";
    #for ($i=0; $i <= $n; $i++) {
    #   if ( !defined $sse[$i] )  {$sse[$i] = $sse[$i-1]};
    #   print "$i $sse[$i]\n"
    #}
    #for ($i=0; $i <= $n; $i++) {
    
    # print "$i $sse[$i]\n"
    #}
    #undef @sse;
    #    $num++;
    
    close(FILE);
}

#print "$header";
#for ($i=0; $i < @tsse; $i++) {
#    $tsse[$i] = $tsse[$i] / $count[$i];
#    print "$i\t$tsse[$i]\n";
#}


#foreach (@ARGV){
#    print;
#}
