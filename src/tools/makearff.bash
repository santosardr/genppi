#!/bin/bash
proteina1=$1
proteina2=$2
db=0
n=1
error_and_exit() {
  echo "$1"  # Output the error message passed as the first argument
  exit 1     # Exit with a status code of 1 to indicate an error
}
progress_bar() {
    porcentagem=$1
    total=100
    # Calcula quantos '#' devem ser impressos
    num_barras=$((porcentagem*total/100))
    num_espacos=$((total-num_barras))
    # Cria a barra de progresso
    barras=$(printf '▓%.0s' $(seq -s ' ' 1 $num_barras))
    espacos=$(printf '░%.0s' $(seq -s ' ' 1 $num_espacos))
    # Imprime a barra de progresso
    printf '|%s%s| %d%%\r' "$barras" "$espacos" "$porcentagem"
}
show_spinner(){
    local -r pid="${1}"; local -r delay='0.75'; local spinstr='\|/-'; local temp
    echo "Executando o comando blastp. Por favor, aguarde..."
    while kill -0 "${pid}" 2>/dev/null; do
	temp="${spinstr#?}"
	printf " [%c]  " "${spinstr}"
	spinstr=${temp}${spinstr%"${temp}"}
	sleep "${delay}"
	printf "\b\b\b\b\b\b"
    done
    echo "Fim do blastp."
    printf "    \b\b\b\b"
}
#validação dos programas e arquivos
if [ $# -lt 2 ]
then
    echo "Parâmetros insuficientes. Exemplos:"; echo "./makearff query.faa subject.faa"; echo "./makearff query.faa subject.faa db"
    exit
fi
if [ $# -eq 3 ]
then
    echo "Assumindo subject como DataBase"; db=1
fi
if [ $# -eq 4 ]
then 
	neg=$4
else
	neg=10	
fi
echo "Porcentagem de negativos $neg";

# Array of protein files to check
files=("$proteina1" "$proteina2" "propensity.dat")
# Loop through the protein files and check if they exist
for file in "${files[@]}"; do
  if [ ! -f "$file" ]; then
    error_and_exit "Arquivo não encontrado: $file"
  fi
done

dir=$(dirname "$proteina1")

commands=("valifasta" "makeblastdb" "blastp" "countfasta" "features")
for cmd in "${commands[@]}"; do
  if [ ! -r "$(which $cmd)" ]; then
    error_and_exit "$cmd not found"
  fi
done
#começo da execução do script
if [ proteina1 -a proteina2 ]
then
    if [ $db -eq 0 ]; then	
	valifasta -i "$proteina1" -o "$proteina1"
	valifasta -i "$proteina2" -o "$proteina2"
	blastp -query "$proteina1" -subject "$proteina2"  -out $dir/blast_protein -evalue 0.000001 -outfmt 6 &
    	show_spinner "$!"
    else
	makeblastdb -dbtype prot -in $proteina2 -hash_index -out $dir/subjectdb
	blastp -query "$proteina1" -db $dir/subjectdb  -out $dir/blast_protein -evalue 0.000001 -outfmt 6 -num_threads 8 &
	show_spinner "$!"
    fi 
    progress_bar 15
    if [ "$proteina1" != "$proteina2" ]
    then
	cat $proteina2 $proteina1 >> $dir/protein.faa
    else
	cat $proteina1 >> $dir/protein.faa
    fi 
    progress_bar 18
    awk '{ if ($1 != $2) { printf("%s\n", $0) } }' $dir/blast_protein > $dir/blast_protein.distintas
    progress_bar 20
    countfasta $dir/protein.faa > $dir/protein.size
    progress_bar 25
    awk -v dir="$dir" '{if ( int($3)>=65 && int ($8-$7)>= int(0.90*($10-$9)) ){ cmd = "grep " $1 " " dir "/protein.size | cut -f 2"; if ( (cmd | getline search_result) > 0 ) { split(search_result, size, "\t"); if (int($8-$7) >= int(0.90*size[1])){ printf("%s\n", $0); } } close(cmd); }}' $dir/blast_protein.distintas > $dir/blast_protein_positivo
    progress_bar 30
    awk -v neg="$neg" 'BEGIN { srand();neg=neg/100}{ rand_val = rand();if (int($3) < 65 && rand_val < neg) { printf("%s\n", $0) } }' $dir/blast_protein.distintas | grep -v ^[0-9] > $dir/blast_protein_negativo
    progress_bar 35
    awk '{ printf("%s\t%s\n", $1, $2)  }' $dir/blast_protein_positivo > $dir/protein_positivo.pares
    progress_bar 40
    awk '{ printf("%s\t%s\n", $1, $2)  }' $dir/blast_protein_negativo > $dir/protein_negativo.pares
    progress_bar 45
    awk '{ key = $1 < $2 ? $1 FS $2 : $2 FS $1 }!sen[key]++' $dir/protein_positivo.pares > $dir/protein_positivo_sem_repeticoes.pares
    progress_bar 50
    awk '{ key = $1 < $2 ? $1 FS $2 : $2 FS $1 }!sen[key]++' $dir/protein_negativo.pares > $dir/protein_negativo_sem_repeticoes.pares
    progress_bar 55
    features $dir/protein.faa > $dir/protein.features
    progress_bar 60
    line=$(head -n 1 $dir/protein.features)
    echo -n -e "\t" > $dir/cabecalho.features
    echo -n $line " " | sed 's/ /1\t/g' | sed 's/1\t$//g' >> $dir/cabecalho.features
    echo -n -e $line "\t" | sed 's/ /2\t/g' | sed 's/\t$/\n/g' >> $dir/cabecalho.features
    progress_bar 65
    cp $dir/cabecalho.features $dir/protein_positivo.features
    progress_bar 70
    cp $dir/cabecalho.features $dir/protein_negativo.features
    progress_bar 75
    declare -A arr2
    while IFS=$'\t' read -r -a linha
    do
	chave=${linha[0]}
	seq=$(IFS=$'\t'; echo "${linha[*]:1}")
	arr2[$chave]=$seq
    done < $dir/protein.features
    progress_bar 80
    while IFS=$'\t' read -r -a arr1
    do
	chave1=${arr1[0]}
	chave2=${arr1[1]}
	seq1=${arr2[$chave1]}
	seq2=${arr2[$chave2]}
	printf "%s - %s\t%s\t%s\n" "$chave1" "$chave2" "$seq1" "$seq2" >> $dir/protein_positivo.features
    done < $dir/protein_positivo_sem_repeticoes.pares
    progress_bar 85
    declare -A arr2
    while IFS=$'\t' read -r -a linha
    do
	chave=${linha[0]}
	seq=$(IFS=$'\t'; echo "${linha[*]:1}")
	arr2[$chave]=$seq
    done < $dir/protein.features
    progress_bar 90 
    while IFS=$'\t' read -r -a arr1
    do
	chave1=${arr1[0]}
	chave2=${arr1[1]}
	seq1=${arr2[$chave1]}
	seq2=${arr2[$chave2]}
	printf "%s - %s\t%s\t%s\n" "$chave1" "$chave2" "$seq1" "$seq2" >> $dir/protein_negativo.features
    done < $dir/protein_negativo_sem_repeticoes.pares
    progress_bar 95
    for file in $dir/blast_protein $dir/blast_protein.distintas $dir/blast_protein_negativo $dir/blast_protein_positivo $dir/cabecalho.features $dir/protein.faa $dir/protein.features $dir/protein_negativo.pares $dir/protein_negativo_sem_repeticoes.pares $dir/protein_positivo.pares $dir/protein_positivo_sem_repeticoes.pares $dir/protein.size $dir/subjectdb.p*; do
	if [ -r $file ];then
            rm $file
	fi
    done
    if [ -r $dir/protein_positivo.features -a -r $dir/protein_negativo.features ];then
	head -n 1 $dir/protein_positivo.features > $dir/weka.attributes
	cp $dir/protein_positivo.features $dir/wekapos.arff
	cp $dir/protein_negativo.features $dir/wekaneg.arff
	sed -i '1d' $dir/weka*.arff
	sed -i "s/\([a-zA-Z0-9]\+\)/@attribute \1 numeric#/g" $dir/weka.attributes
	tr '#' '\n' < $dir/weka.attributes > $dir/weka.attributes2
	tr -d '\t' < $dir/weka.attributes2 > $dir/weka.attributes
	cut -f 2- $dir/wekapos.arff > $dir/wekapos.arff2
	cut -f 2- $dir/wekaneg.arff > $dir/wekaneg.arff2
	sed -i "s/\t/,/g" $dir/weka*.arff2
	sed -i "s/$/,POSITIVE/g" $dir/wekapos.arff2
	sed -i "s/$/,NEGATIVE/g" $dir/wekaneg.arff2
	echo '@relation similar' > $dir/similar.arff
	cat $dir/weka.attributes  >> $dir/similar.arff
	echo '@attribute class {POSITIVE,NEGATIVE}' >> $dir/similar.arff
	echo '@data' >> $dir/similar.arff
	progress_bar 97
	cat $dir/similar.arff > $dir/similar_sem_wekaneg-wekapos.arff
	progress_bar 98
	cat $dir/wekaneg.arff2 $dir/wekapos.arff2  >> $dir/similar.arff
	progress_bar 99
	for file in $dir/weka.attributes $dir/wekapos.arff $dir/wekaneg.arff $dir/weka.attributes2; do
	    if [ -r $file ];then
		rm $file
	    fi
	done
	mv $dir/wekapos.arff2 $dir/wekapos.arff
	mv $dir/wekaneg.arff2 $dir/wekaneg.arff
	progress_bar 100
	echo -e "\nPronto, aqui está seus arquivos desejados, o arquivo final se chama similar.arff."
    fi
else
    echo "Para continuar voce precisa ter todos esses arquivos no seu diretório: "
fi
