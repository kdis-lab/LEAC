#run:
# awk -f movement_libras_trans.awk -F ',' -v OFS=',' movement_libras.data  > movement_libras_trasn.data
{
    for(i=1;i< NF;i++){ t1+=$i; i++; t2+=$i };
    xd = 0.5 - 2.0 * t1/(NF-1);
    yd = 0.5 -2.0 * t2/(NF-1);
    t1=0; t2=0
    for(i=1;i< NF;i++){ $i +=  xd;  i++; $i += yd };
    for(i=1; i<=NF; i++) printf "%s",$i (i==NF?ORS:OFS)     
}
