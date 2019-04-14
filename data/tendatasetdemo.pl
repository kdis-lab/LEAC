#!/usr/bin/perl
# /*! \file tendatasetdemo.pl
#  *
#  * \brief  Crea los archivos de 10 dataset para friedman test with k-fold 
#  *
#  * \details  This file is part of the LEAC.\n\n
#  * \version 1.0
#  * \date 2015-2017
#  * \authors Hermes Robles-Berumen <hermes@uaz.edu.mx>\n Sebastian Ventura <sventura@uco.es>\n Amelia Zafra <azafra@uco.es>\n <a href="http://www.uco.es/kdis/">KDIS</a>
#  * \copyright <a href="https://www.gnu.org/licenses/gpl-3.0.en.html">GPLv3</a> license
#  */
#
#
#  Use:
#       perl tendatasetdemo.pl
#

use strict;
use warnings;

use LWP::Simple;
use File::Copy;

sub insertheader {
    # Get passed arguments
    my ($namefile, $header) = @_;

    open(FILE,$namefile) || die "can't open file for read\n"; 
    my @lines=<FILE>;
    close(FILE);
    open(FILE,">$namefile")|| die "can't open file for write\n";
    print FILE "$header\n" ;
    foreach my $line (@lines) {
	print FILE $line;
    }#end foreach
    close(FILE);
}

#IRIS DATASET
my $url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data';
my $file = 'iris.data';
my $header = 'SepalLength,SepalWidth,PetalLength,PetalWidth,Class';
getstore($url,$file);
insertheader($file,$header);
print "Download $file OK\n";
system('../eac/stdvar_milligan_cooper1988 -i iris.data -a "1-4" -c 5 -h yes --std-var Z5 > iris_z5.data');
print "Stdvar $file OK\n";
system("Rscript kfold.R iris_z5.data");
print "k-fold $file OK\n";

#WINE DATASET
$url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/wine/wine.data';
$file = 'wine.data';
$header = 'Alcohol,MalicAcid,Ash,AlcalinityAsh,Magnesium,TotalPhenols,Flavanoids,NonflavanoidPhenols,Proanthocyanins,ColorIntensity,Hue,OD280_OD315,Proline';
getstore($url, $file);
insertheader($file,$header);
print "Download $file OK\n";
system('../eac/stdvar_milligan_cooper1988 -i wine.data -a "2-13" -c 1 -h yes --std-var Z5 > wine_z5.data');
print "Stdvar $file OK\n";
system("Rscript kfold.R wine_z5.data");
print "k-fold $file OK\n";


#SONAR DATASET
$url = 'http://archive.ics.uci.edu/ml/machine-learning-databases/undocumented/connectionist-bench/sonar/sonar.all-data';
$file = 'sonar.all-data';
$header = 'X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,X21,X22,X23,X24,X25,X26,X27,X28,X29,X30,X31,X32,X33,X34,X35,X36,X37,X38,X39,X40,X41,X42,X43,X44,X45,X46,X47,X48,X49,X50,X51,X52,X53,X54,X55,X56,X57,X58,X59,X60,Class';
getstore($url, $file);
insertheader($file,$header);
print "Download $file OK\n";
system('../eac/stdvar_milligan_cooper1988 -i sonar.all-data -a "1-60" -c 61 -h yes --std-var Z5 > sonar_z5.data');
print "Stdvar $file OK\n";
system("Rscript kfold.R sonar_z5.data");
print "k-fold $file OK\n";


#GLASS DATASET 4
$url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/glass/glass.data';
$file = 'glass.data';
$header = 'Id,RI,Na,Mg,Al,Si,K,Ca,Ba,Fe,TypeGlass';
getstore($url, $file);
insertheader($file,$header);
print "Download $file OK\n";
system('../eac/stdvar_milligan_cooper1988 -i glass.data -a "2-10" -c 11 -h yes --std-var Z5 > glass_z5.data');
print "Stdvar $file OK\n";
system("Rscript kfold.R glass_z5.data");
print "k-fold $file OK\n";


#SPECTF DATASET 5
$header = 'DIAGNOSIS,F1,F2,F3,F4,F5,F6,F7,F8,F9,F10,F11,F12,F13,F14,F15,F16,F17,F18,F19,F20,F21,F22';
$url  = 'https://archive.ics.uci.edu/ml/machine-learning-databases/spect/SPECT.train';
$file = 'SPECT.train';
getstore($url, $file);
$url  = 'https://archive.ics.uci.edu/ml/machine-learning-databases/spect/SPECT.test';
$file = 'SPECT.test';
getstore($url, $file);
copy "SPECT.train", "spect.data"; # or die "Copy failed: $!";
open FILEA, '>> spect.data' or die $!;
open FILEB, '< SPECT.test' or die $!;
print FILEA <FILEB>;
insertheader('spect.data',$header);
print "Download $file OK\n";
system('../eac/stdvar_milligan_cooper1988 -i spect.data -a "2-23" -c 1 -h yes --std-var Z5 > spect_z5.data');
print "Stdvar $file OK\n";
system("Rscript kfold.R spect_z5.data");
print "k-fold $file OK\n";

#HABERMAN DATASET 6
$url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/haberman/haberman.data';
$file = 'haberman.data';
$header = 'Age,PatientsYear,Nodes,Survival';
getstore($url, $file);
insertheader($file,$header);
print "Download $file OK\n";
system('../eac/stdvar_milligan_cooper1988 -i haberman.data -a "1-3" -c 4 -h yes --std-var Z5 > haberman_z5.data');
print "Stdvar $file OK\n";
system("Rscript kfold.R haberman_z5.data");
print "k-fold $file OK\n";

#ECOLI DATASET 7
$url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/ecoli/ecoli.data';
$file = 'ecoli.data';
$header = 'name mcg gvh lip chg aac alm1 alm2 class';
getstore($url, $file);
insertheader($file,$header);
print "Download $file OK\n";
system('../eac/stdvar_milligan_cooper1988 -i ecoli.data -d " " -a "2-8" -c 9 -h yes --std-var Z5 > ecoli_z5.data');
print "Stdvar $file OK\n";
system("Rscript kfold.R ecoli_z5.data");
print "k-fold $file OK\n";

#IONOSPHERE DATASET 8
$url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/ionosphere/ionosphere.data';
$file = 'ionosphere.data';
$header = 'X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20,X21,X22,X23,X24,X25,X26,X27,X28,X29,X30,X31,X32,X33,X34,Class';
getstore($url, $file);
insertheader($file,$header);
print "Download $file OK\n";
system('../eac/stdvar_milligan_cooper1988 -i ionosphere.data -a "1,3-34" -c 35 -h yes --std-var Z5 > ionosphere_z5.data');
print "Stdvar $file OK\n";
system("Rscript kfold.R ionosphere_z5.data");
print "k-fold $file OK\n";

#HAYES-ROTH DATASET 9
$header = 'name,hobby,age,educational,marital,class';
$url  = 'https://archive.ics.uci.edu/ml/machine-learning-databases/hayes-roth/hayes-roth.data';
$file = 'hayes-roth.data';
getstore($url, $file);
print "Download $file OK\n";
$url  = 'https://archive.ics.uci.edu/ml/machine-learning-databases/hayes-roth/hayes-roth.test';
$file = 'hayes-roth.test';
getstore($url, $file);
#ADD A COLUMN TO hayes-roth.test
open(FILE,'hayes-roth.test') || die "can't open file for read\n"; 
my @lines=<FILE>;
close(FILE);
open(FILE,">hayes-roth.test")|| die "can't open file for write\n";
foreach my $line (@lines) {
  print FILE "-1,$line";
}#end foreach
close(FILE);
#
copy "hayes-roth.data", "hayes_roth.data"; # or die "Copy failed: $!";
open FILEA, '>> hayes_roth.data' or die $!;
open FILEB, '< hayes-roth.test' or die $!;
print FILEA <FILEB>;
insertheader('hayes_roth.data',$header);
print "Download $file OK\n";
system('../eac/stdvar_milligan_cooper1988 -i hayes_roth.data -a "3-5" -c 6 -h yes --std-var Z5 > hayes_roth_z5.data');
print "Stdvar $file OK\n";
system("Rscript kfold.R hayes_roth_z5.data");
print "k-fold $file OK\n";


#ZOO DATASET 10
$url = 'https://archive.ics.uci.edu/ml/machine-learning-databases/zoo/zoo.data';
$file = 'zoo.data';
$header = 'animalname,hair,feathers,eggs,milk,airborne,aquatic,predator,toothed,backbone,breathes,venomous,fins,legs,tail,domestic,catsize,type';
getstore($url, $file);
insertheader($file,$header);
print "Download $file OK\n";
system('../eac/stdvar_milligan_cooper1988 -i zoo.data -a "2-17" -c 18 -h yes --std-var Z5 > zoo_z5.data');
print "Stdvar $file OK\n";
system("Rscript kfold.R zoo_z5.data");
print "k-fold $file OK\n";


