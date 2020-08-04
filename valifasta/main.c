#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "validador.h"


int main(int argc, char *argv[ ])
{
    char modo_leitura[1] = "r";
    char modo_escrita_continua[1] = "a";
    char modo_escrita[1] = "w";
    char caminho_arquivo_leitura[300];
    char caminho_arquivo_escrita[300];
    char caminho_arquivo_log[300];
    char comandoInsert[3] = "-i";
    char comandoOut[3] = "-o";
    char comandoLog[3] = "-l";
    char comandoHelp[6] = "-help";
    char comandoHELP[6] = "-HELP";

    if(argc < 2){
        printf("Missing arguments: -i (for insert) name_of_file (to be processeced)!\n");
        getchar();
        return 1;
    }

    int i, comandos[10] = {0,0,0,0,0,0,0,0,0,0};
    char ultlido[100] = " ";

    for(i=1;i<argc;i++)
    {
        if(!strcmp(argv[i], comandoInsert))
        {
            comandos[0]++;
        }
        else
        if(!strcmp(argv[i], comandoOut))
        {
            comandos[1]++;
            geraArqSaida = gerarNomeSaida = 1;
        }
        else
        if(!strcmp(argv[i], comandoLog))
        {
            comandos[2]++;
            geralog = gerarNomeLog = 1;
        }
        else
        if(!strcmp(argv[i], comandoHELP) || !strcmp(argv[i], comandoHelp))
        {
            comandos[3]++;
            printf("Usage: ./main -i [file_for_be_reading] -o [file_of_out] --Result in a file passed\n");
            printf("Usage: ./main -i [file_for_be_reading] -o --Result in a new file (file_of_input_saida.fasta)\n");
            printf("Usage: ./main -i [file_for_be_reading] -o [file_of_out] -l [file_of_log] --Result in a file passed and creating log in a file passed\n");
            printf("Usage: ./main -i [file_for_be_reading] -o -l --Result in a new file (file_of_input_saida.fasta) and creating a file of log (file_of_input_logger.log)\n");
            printf("Usage: ./main -i [file_for_be_reading] --Result on window\n");
            getchar();
            return 1;
        }
        else
        {
            if(!strcmp(ultlido, comandoInsert))
            {
                strcpy(caminho_arquivo_leitura,argv[i]);
            }
            else
            if(!strcmp(ultlido, comandoOut))
            {
                strcpy(caminho_arquivo_escrita, argv[i]);
                gerarNomeSaida = 0;
            }
            else
            if(!strcmp(ultlido, comandoLog))
            {
                strcpy(caminho_arquivo_log,argv[i]);
                gerarNomeLog = 0;
            }
            else{
                printf("There is no definition for the command %s!\n", argv[i]);
                getchar();
                return 1;
            }
        }
        strcpy(ultlido, argv[i]);
    }

    if(!comandos[0] || argc < 3)
    {
        printf("Need a file for input: -i [arquivo_entrada] [opicional]-o [arquivo_saida]\n");
        getchar();
        return 1;
    }

    char* temp;
    int j=0;
    i=0;

    if(gerarNomeSaida){
        ///Alocando memoria para o arquivo de log
        strcpy(caminho_arquivo_escrita, "");
        strcpy(caminho_arquivo_escrita, caminho_arquivo_leitura);
        strcat(caminho_arquivo_escrita, "_saida.fasta");
    }

    if(gerarNomeLog){
        ///Alocando memoria para o arquivo de log
        strcpy(caminho_arquivo_log, "");
        strcpy(caminho_arquivo_log, caminho_arquivo_leitura);
        strcat(caminho_arquivo_log, "_loger.log");
    }

    FILE* log;

    if(geralog){
        if((log = fopen(caminho_arquivo_log, modo_escrita)) == NULL)
        {
            printf("It is not possible to open the file of logs\n");
            return 1;
        }
        fprintf(log, "-------Starting\n\nFile of log is open!\n");
        fclose(log);
    }


    Fila* F;
    cabecalho* C;
    sequencia* S;

    ///Cria a fila para armazenar os dados do arquivo a ser lido
    F = Cria();

    FILE* file;

    char c;
    char* token;


    ///Lê o arquivo e o armazena em file
    if((file = fopen(caminho_arquivo_leitura, modo_leitura)) == NULL)
    {
        if(!geralog)
            return 1;

        log = fopen(caminho_arquivo_log, modo_escrita_continua);
        fprintf(log, "It is not possible to open the file of input! File: %s\n", caminho_arquivo_leitura);
        fprintf(log, "ERROR!\n");
        fclose(log);
        return 1;
    }


    int sizeOfsequence = 10000;
    int numero_proteinas=0, ehcabecalho = 1;

    i=0;
    int x = 0;
    long int linha = 0;
    int erro = 0;
    temp = (char*)malloc(sizeOfsequence*sizeof(char));

    if(geralog){
        log = fopen(caminho_arquivo_log, modo_escrita_continua);

        fprintf(log, "File of input is open! File: %s\n", caminho_arquivo_leitura);
        fclose(log);
    }

    int comecei=1;

    ///Percorre todo o arquivo atribuindo o caracter atual para a variável c
    while( x>=0 )
    {
        x = fgetc(file);
        c = (char)x;

        if(erro)
        {
            while(c != '>')
            {
                x = fgetc(file);
                c = (char)x;
                if(x<0)
                {
                    break;
                }
            }
            erro = 0;
        }

        ///Valida se o primeiro caracter a se ler é '<' ou quebra linha até achar um caracter
        if(comecei)
        {
            comecei--;
            if(c != '>')
            {
                if(geralog){
                    log = fopen(caminho_arquivo_log, modo_escrita_continua);
                    fprintf(log, "Não há definição de cabeçalho, é necessário que cada cabeçalho inicie com '>'");
                    fclose(log);
                }
                return 1;
            }
        }

        ///Coloca o caracter em letra maiúscula se for uma letra e se estiver na sequência
        if(x >= 97 && x <= 122 && !ehcabecalho)
        {
            c-=32;
            x-=32;
        }

        if(c == '>')
            ehcabecalho = 1;

        if(linha && c == '>' || x < 0)
        {
            temp[i] = '\0';
            ehcabecalho = 1;
            S = (sequencia*)malloc(sizeof(sequencia));
            S->sequencia = (char*)malloc((strlen(temp)+1) * sizeof(char));
            strcpy(S->sequencia, temp);
            S->tamanho = strlen(S->sequencia);

            Entra(F, C, S);
            free(temp);
            temp = (char*)malloc(sizeOfsequence*sizeof(char));
            i=0;
        }

        if(c == '\n')
        {
            linha++;

            if(ehcabecalho)
            {
                temp[i] = '\0';
                ehcabecalho = 0;

                C = (cabecalho*)malloc(sizeof(cabecalho*));
                C->nome = (char*)malloc((strlen(temp)+1) * sizeof(char));
		if(strlen(temp)<64)
		  strcpy(C->nome, temp);
		else
		  strncpy(C->nome, temp,64);
                free(temp);
                temp = (char*)malloc(sizeOfsequence*sizeof(char));
                i=0;
            }
        }
        else
        {
            if(i < (sizeOfsequence-1))
            {
                if(!ehcabecalho && x >= 65 && x <= 90)
                {
                    temp[i] = c;
                    i++;
                }
                else
                {
                    if(ehcabecalho)
                    {
                        temp[i] = c;
                        i++;
                    }
                }
            }
            else
            {
                if(geralog)
                {
                    log = fopen(caminho_arquivo_log, modo_escrita_continua);
                    fprintf(log, "Sequencia de aminoácidos com mais de %d aminoacidos, impossível de processar\n", sizeOfsequence);
                    fclose(log);
                }
                erro = 1;
            }
        }
    }

    SetaValores(caminho_arquivo_escrita, caminho_arquivo_leitura, caminho_arquivo_log);
    if(geralog)
    {
        log = fopen(caminho_arquivo_log, modo_escrita_continua);
        fprintf(log, "Comecar validacoes\n");
        fclose(log);
    }


    /// Valida cabeçalho dos aminoácidos da Fila
    int retorno = ValidaCabecalhoAminoacido(F, log);
    if(!retorno)
    {
        if(geralog)
        {
            log = fopen(caminho_arquivo_log, modo_escrita_continua);
            fprintf(log, "Não há cabeçalhos para descarregar\n");
            fprintf(log, "ERRO!\n");
            fclose(log);
        }
        return 1;
    }
    retorno = 0;

    if(geralog)
    {
        log = fopen(caminho_arquivo_log, modo_escrita_continua);
        fprintf(log, "Cabecalhos Validados\n\n");
        fclose(log);
    }


    /// Valida sequência dos aminoácidos da Fila
    retorno = ValidaSequenciaAminoacido(F, log);

    if(!retorno)
    {
        if(geralog)
        {
            log = fopen(caminho_arquivo_log, modo_escrita_continua);
            fprintf(log, "Não há sequencias para validar\n");
            fprintf(log, "ERRO!\n");
            fclose(log);
        }
        return 1;
    }

    if(geralog){
        log = fopen(caminho_arquivo_log, modo_escrita_continua);
        fprintf(log, "Sequencias Validadas\n");
        fclose(log);
    }
    if (file) fclose(file);

    if(geraArqSaida)
    {
        if((file = fopen(caminho_arquivo_escrita, modo_escrita)) == NULL)
        {
            if(geralog)
            {
                log = fopen(caminho_arquivo_log, modo_escrita_continua);
                fprintf(log, "Não foi possível abrir o arquivo de saida! Aquivo: %s\n", caminho_arquivo_escrita);
                fprintf(log, "ERRO!\n");
                fclose(log);
            }
            return 1;
        }
        //log = fopen(caminho_arquivo_log, modo_escrita_continua);
        /// Descarrega no arquivo de saída todo o arquivo
        if(DescarregaFilaEmArquivo(F, file, log))
            return 1;
        //fclose(log);
        if (file) fclose(file);
    }
    else
    {
        if(DescarregaFilaEmTela(F, log))
            return 1;
    }


    if(geralog)
    {
        log = fopen(caminho_arquivo_log, modo_escrita_continua);

        fprintf(log, "File: %s finished. \nLog: %s Output: %s\n\n\n", caminho_arquivo_leitura, caminho_arquivo_log, caminho_arquivo_escrita);
        fprintf(log, "Finished!");

        fclose(log);
    }


    Destroi(F);

    return 0;
}
