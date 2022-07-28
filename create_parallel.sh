#!/bin/bash

trait=$1
folder_name=$2

for i in {1..64}
do 
  echo "#!/bin/bash" > "/nethome/kcni/jxia/chromatin-accessibility/${folder_name}/${trait}_${i}.sh"
  echo "Rscript /nethome/kcni/jxia/chromatin-accessibility/parallel.R ${i}" >> "/nethome/kcni/jxia/chromatin-accessibility/${folder_name}/${trait}_${i}.sh"
  echo "/nethome/kcni/jxia/chromatin-accessibility/${folder_name}/${trait}_${i}.sh" >> "/nethome/kcni/jxia/chromatin-accessibility/${folder_name}/script_list.txt" 
  chmod +x /nethome/kcni/jxia/chromatin-accessibility/${folder_name}/${trait}_${i}.sh
  echo "File ${trait}_${i}.sh generated"  
done
