#run:
# awk -f zoo_binary.awk -F ',' -v OFS=','  zoo.data > zoo_bin.csv
BEGIN {

    h1 = "animalname";
    h2 = "hair";
    h3 = "feathers";
    h4 = "eggs";
    h5 = "milk";
    h6 = "airborne";
    h7 = "aquatic";
    h8 = "predator";
    h9 = "toothed";
    h10 = "backbone";
    h11 = "breathes";
    h12 = "venomous";
    h13 = "fins";
    h14 = "legs_0,legs_2,legs_4,legs_5,legs_6,legs_8";
    h15 = "tail";
    h16 = "domestic";
    h17 = "catsize";
    h18 = "type";
  
    print h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h14,h15,h16,h17,h18;
}

{
# legs:Numeric (set of values: {0,2,4,5,6,8}) 
    if ( $14 == 0)
	$14 = "1,0,0,0,0,0";
    else if ($14 == 2)
	$14 = "0,1,0,0,0,0";
    else if ($14 == 4)
	$14 = "0,0,1,0,0,0";
    else if ($14 == 5)
	$14 = "0,0,0,1,0,0";
    else if ($14 == 6)
	$14 = "0,0,0,0,1,0";
    else if ( $14 == 8)
	$14 = "0,0,0,0,0,1";
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18
}
