#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SEQUENCE_NAME 1000
int main (int argc, char *argv[] )
{
	FILE * pFile, * pOut = NULL;
	char *sequence_name, c;
	printf("Separando sequencias fasta do arquivo %s:\n", argv[1]);
	sequence_name = (char*) malloc(MAX_SEQUENCE_NAME);
	*sequence_name='\0';
	//Abre o arquivo de origem dos dados
	pFile = (FILE * ) fopen (argv[1],"r");
	//Se o arquivo eh valido proceda ...
	if (pFile!=NULL){
		do {//leia todos os caracteres do arquivo um a por um
			c = getc (pFile);
			if (c == '>'){//se encontrar um caracter de início de sequencia fasta 
				//fecha algum arquivo de saida que pode estar aberto
				if (pOut!=NULL) fclose(pOut);
				//inicia a variavel com o nome da sequencia
				sequence_name[0]='\0';
				//declara e inicia variaveis que dizem o tamanho do nome e
				//onde estah o primeiro caractere separador de texto dentro
				//do nome da sequencia. Nesse ponto ignoramos o resto do texto
				//para o nome de nossa copia da sequencia
				int sn_size, sn_stop;
				sn_size=0; sn_stop=0;
	 			//Le o nome da sequencia fasta caracter por caracter
				do{
					c = getc (pFile);
					sequence_name[sn_size]=c;
					sn_size=sn_size+1;
					//Se encontrar um separador de texto marca onde foi
					if( isspace(c) && !sn_stop ) 
						sn_stop = sn_size;
						//continue a leitura ateh o fim da linha
				}while (c != EOF && c != '\n' && c != '\r');
				//Se encontrou um caracter separador termine o nome da sequencia nesse ponto
				if (sn_stop>0) sn_size = sn_stop;
					sequence_name[sn_size-1]='\0';
					printf("%s ",sequence_name);
					//cria outro arquivo para armazenar a nova sequencia
					pOut = fopen( sequence_name, "w");
					if (pOut!=NULL){
						//Escreve o simbolo de inicio de seq fasta e nome da seq
						fputc('>', pOut);
						fputs(sequence_name, pOut);
					}//if pOut!=NULL
				}//if c == '>'
			//Apos verificar se era um caractere de inicio de sequencia fasta
			//acrescente um caracter lido (diferente de '>') ao arquivo da sequencia
			if (pOut!=NULL){
				if (c != '>' && isascii(c) ) fputc(c, pOut);
			}//if pOut!=NULL
		} while (c != EOF);//fim do loop principal
		fclose (pFile);//fecha o arquivo de origem
	}//if (pFile!=NULL)
	free(sequence_name);//libera memoria
	return 0;
}//fim programa
