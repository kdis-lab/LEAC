#run:
# awk -f movement_libras_trans.awk -F "," -v OFS="," movement_libras.data  > movement_libras_trans.data
{
    for(i=1;i< NF;i++){ sumcoord1+=$i; i++; sumcoord2+=$i };
    numpoint = (NF-1) / 2.0;
    xd = 0.5 - sumcoord1/numpoint;
    yd = 0.5 - sumcoord2/numpoint;
    sumcoord1=0; sumcoord2=0;
    for(i=1;i< NF;i++){ $i +=  xd;  i++; $i += yd };
    for(i=1; i<=NF; i++) printf "%s",$i (i==NF?ORS:OFS)     
}
