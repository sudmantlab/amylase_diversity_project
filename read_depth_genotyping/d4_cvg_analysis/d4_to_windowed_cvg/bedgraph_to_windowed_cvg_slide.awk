#takes in bedgraph
#then outputs sliding windows
#usage cat {bedgraph} | awk -v w=1000 -v s=200 -f graph_to_windowed_cvg_slide.awk 
BEGIN{
    idx = 1; 
    #w = 1000; ##INCLUDE AS A variabler
    #slide = 1000; 
    sum = 0;
    for (j=0; j<w; j++){
        arr[j] = 0;
    }
    #printf("#w:%d s:%d",w,slide)
}
#FOR EACH LINE
{
    
    contig=$1 
    cvg=$4
    sample=$5
    for (i = $2; i<$3; i++){
        mod = idx%w;
        sum = sum-arr[mod]
        sum = sum+cvg;
        arr[mod] = cvg

        start = i;
        end = i+1;

        if (((idx % s)==0)&&(idx>=w)){
            printf("%s\t%d\t%d\t%f\t%s\n",contig,starts_arr[mod],end,sum/w,sample);
        }
        starts_arr[mod] = end 
        idx = idx+1;
    }
}
