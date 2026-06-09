input_file="./E.-coli-RBS-Prediction-and-Optimization/50nt_upstream_30nt_ORF.fasta"
output_dir="./50nt_2nd"
if [ ! -d "$output_dir" ]; then
mkdir "$output_dir"
fi
while read -r line1 && read -r line2; do
gene_id=$(echo "$line1" | cut -d "|" -f1-2 | tr -d '>')
dot_bracket=$(echo "$line2" | RNAfold --noPS --noLP)
echo "$dot_bracket" > "$output_dir/$gene_id.fold"
done < "$input_file"
