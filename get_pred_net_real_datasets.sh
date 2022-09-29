
#inputs_dir="inputs/datasets"
results_dir="/home/mengxu/scGRN-L0/output/"
postfix="-string-700"

# run this from your output datasets dir 
#out_dir="/data/inputs/single-cell/RNMethods/"
out_dir="/home/mengxu/scGRN-L0/output/"
#num_edges=100; 
#base_file="rankedEdges.csv"; 
base_file="output_net_hESC_L0.csv"; 
for f in `find $results_dir -maxdepth 4 -name "$base_file"`; do 
    echo $f
    base_dir=`dirname $f`
    #inputs_dir=`dirname $(dirname $f) | sed "s/output/inputs/g"`
    inputs_dir="output/"
    #echo $inputs_dir
    # get the number of edges from the original string network. The header line is ok because the rankedEdges files also have a header
    num_edges=`cat ${inputs_dir}/refNetwork${postfix}.csv | wc -l`
    #echo $num_edges

    mkdir -p "${out_dir}/${base_dir}"
    out_file="${out_dir}/${base_dir}/rankedEdges${postfix}.csv"
    tmp_file="${out_dir}/${base_dir}/rankedEdges${postfix}-noself.csv"
    echo "Writing subnetwork to: ${out_file}"

    # make sure there are no self edges
    cat $f | awk 'BEGIN{FS=OFS="\t";} {if ($1 != $2) {print $0}}' > $tmp_file
    
    score=`head -n $num_edges $tmp_file | tail -n 1 | cut -f 3`; 
    #echo $score

    # now get the subnetwork
    cat $tmp_file | awk -v score="$score" 'function abs(v) {return v < 0 ? -v : v} BEGIN{FS=OFS="\t";} {if (abs($3) >= abs(score)) {print $0}}' > $out_file
    wc -l $out_file
    break
done
