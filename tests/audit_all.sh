#!/bin/bash
cd ~/projects/hit_validation
mkdir -p audit_outputs

ls -d 04_data/campaigns/*/ | while read d; do
    [ -n "$(ls $d/ligands/*.sdf 2>/dev/null)" ] && echo "$d"
done > /tmp/audit_dirs.txt

cat /tmp/audit_dirs.txt | xargs -P 8 -I DIR bash -c '
d="DIR"
camp=$(basename "$d")
out="audit_outputs/audit_${camp}.csv"
python audit_protonation_states_v3.py "$d/ligands" --ph 6.3 --out "$out" > /dev/null 2>&1
if [ -f "$out" ]; then
    impar=$(python -c "import pandas as pd; df=pd.read_csv(\"$out\"); print((df[\"parity_ok\"]==False).sum() if \"parity_ok\" in df.columns else \"?\")")
    total=$(($(wc -l < "$out") - 1))
    printf "  %-50s  total=%-3s  impar=%s\n" "$camp" "$total" "$impar"
else
    echo "  $camp  FALLÓ"
fi
'