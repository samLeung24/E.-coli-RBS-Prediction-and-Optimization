
for file in individual_gff/*.gff; do
output_file="individual_csv/$(basename "$file" .gff).csv"
awk -F '\t' '$3=="CDS"{print $9}' "$file" > "$output_file"
done

for file in individual_gff/*.gff; do
output_file="individual_csv/$(basename "$file" .gff).csv"
awk -F '\t' '$3=="CDS"{print $1","$4","$5","$7","$9}' "$file" > "$output_file"
done