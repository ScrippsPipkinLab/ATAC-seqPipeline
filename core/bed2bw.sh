module load ucsc_tools;
for file in /gpfs/home/snagaraja/ATACseqPipeline/data/macs2/*.bdg
do  
    base=`basename $file .bdg`;
    newbdg=`basename $file _clipped.bdg`;
    echo $base;
    bedClip $file /gpfs/home/snagaraja/ATACseqPipeline/core/GCF_000001635.27_GRCm39_genomic.size $base.new.bed;
    # echo "Done";
    bedGraphToBigWig $newbdg /gpfs/home/snagaraja/ATACseqPipeline/core/GCF_000001635.27_GRCm39_genomic.size $base.bw;
done
# for file in *.bdg
# do 
#     base=`basename $file .bdg`;
#     echo $base.bw;
# done

# bedGraphToBigWig 