#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "validador.h"

char modo_leitura[1] = "r";
char modo_escrita_continua[1] = "a";
char modo_escrita[1] = "w";
char caminho_arquivo_leitura[300];
char caminho_arquivo_escrita[300];
char caminho_arquivo_log[300];

void Renomear(Fila*);

struct fila
{
    no *inicio, *fim;
    int total;
};

Fila* Cria()
{
    Fila *F = (Fila*)malloc(sizeof(Fila));
    F->inicio = NULL;
    F->fim = NULL;
    F->total=0;
    return F;
}

void Esvazia(Fila *F)
{
    no *ndel, *nextno;
    nextno = F->inicio;

    while (nextno !=NULL)
    {
        ndel = nextno;
        if(ndel->Cabecalho->aminoacidos != NULL)
        {
            free(ndel->Cabecalho->aminoacidos);
        }
        if(ndel->Cabecalho != NULL)
        {
            free(ndel->Cabecalho);
        }
        nextno = nextno->prox;
        free(ndel);
    }

    F->inicio=NULL;
    F->fim=NULL;
    F->total=0;
}


void Destroi(Fila *F)
{
    Esvazia(F);
    free(F);
}

void Entra(Fila *F, cabecalho *Cab, sequencia *Seq)
{
    no *p =(no*) malloc(sizeof(no));
    F->total++;
    Cab->aminoacidos = Seq;
    p->Cabecalho= Cab;
    p->prox=NULL;
    p->correto=1;

    if (F->inicio==NULL)
        F->inicio=p;
    else
        F->fim->prox=p;

    F->fim=p;
}



void Sai(Fila *F, cabecalho *Cab, sequencia *Seq)
{
    no *p;
    if (F->total!=0)
    {
        F->total--;
        p=F->inicio;
        Cab=p->Cabecalho;
        Seq=p->Cabecalho->aminoacidos;
        F->inicio=p->prox;
        if (F->inicio == NULL)
            F->fim = NULL;
        free(p);
    }
}

void Inicio(Fila* F, cabecalho *Cab)
{
    if (!Vazia(F))
    {
        Cab=F->inicio->Cabecalho;
    }
}

int Vazia(Fila *F)
{
    return F->total==0;
}

int Tamanho(Fila *F)
{
    return F->total;
}

int DescarregaFilaEmArquivo(Fila* F, FILE* arquivo, FILE* log)
{
  if(arquivo == NULL)
    {
        return 1;
    }

    no* capsula = F->inicio;

    while(capsula != NULL)
    {
        if(capsula->correto)
        {
            cabecalho* c = capsula->Cabecalho;
            sequencia* s = c->aminoacidos;
            fprintf(arquivo, "%s\n%s\n", c->nome, s->sequencia);
        }
        else
        {
            if(geralog){
                log = fopen(caminho_arquivo_log, modo_escrita_continua);
                fprintf(log, "Sequência do cabeçalho %s está incorreta. ERRO!\n", capsula->Cabecalho->nome);
                fclose(log);
            }
        }
        capsula = capsula->prox;
    }
    return 0;
}

int DescarregaFilaEmTela(Fila* F, FILE* log)
{
    no* capsula = F->inicio;
    while(capsula != NULL)
    {
        if(capsula->correto)
        {
            cabecalho* c = capsula->Cabecalho;
            sequencia* s = c->aminoacidos;
            printf("%s\n%s\n", c->nome, s->sequencia);
        }
        else
        {
            if(geralog)
            {
                log = fopen(caminho_arquivo_log, modo_escrita_continua);
                fprintf(log, "Sequência do cabeçalho %s está incorreta. ERRO!\n", capsula->Cabecalho->nome);
                fclose(log);
            }
        }
        capsula = capsula->prox;
    }
    return 0;
}

void SetaValores(char* c, char* a, char* b)
{
    strcpy(caminho_arquivo_escrita, c);
    strcpy(caminho_arquivo_leitura, a);
    strcpy(caminho_arquivo_log, b);
}

int ValidaUnicidadeCabecalho(no* capsula, char* s, int num, FILE* log)
{
    capsula = capsula->prox;
    char* newcab;
    char c;
    int i = 0, n = num, a=0;

    while(capsula != NULL)
    {
        n=num;
        a=0;

        newcab = (char*)malloc((strlen(capsula->Cabecalho->nome) + 1) * sizeof(char));
        while(n && capsula->Cabecalho->nome[i] != '\0')
        {
            c = capsula->Cabecalho->nome[i];

            if(c == ' ' || c == '|' || c == ',' || c == ';')
            {
                n--;
                newcab[a++] = '_';
                i++;
                continue;
            }

            newcab[a++] = c;
            i++;
        }
        newcab[a] = '\0';
        i=0;

        if(!strcmp(newcab, s))
        {
            return 1;
        }

        capsula = capsula->prox;
    }

    return 0;
}

int ValidaCabecalhoAminoacido(Fila* F, FILE* log)
{
    if(geralog)
    {
        log = fopen(caminho_arquivo_log, modo_escrita_continua);
        fprintf(log, "Iniciando validação dos cabeçalhos da fila.\n\n");
        fclose(log);
    }

    int i=0, erro = 0, n, j=0, a=0, maior =1;
    no* capsula = F->inicio;
    char* frase;
    if(capsula == NULL) return 0;
    while(capsula != NULL)
    {
        i=1;
        erro=a=0;
        while(i < 10)
        {
            n=i;
            j=0;
            a=0;
            frase = (char*)malloc((strlen(capsula->Cabecalho->nome)+1) * sizeof(char));
            while(n && capsula->Cabecalho->nome[j] != '\0')
            {
                if(capsula->Cabecalho->nome[j] == ' ' ||
		   capsula->Cabecalho->nome[j] == '|' ||
		   capsula->Cabecalho->nome[j] == ',' ||
		   capsula->Cabecalho->nome[j] == '.')
                {
                    frase[a++] = '_';
                    n--;
                    j++;
                    continue;
                }
		  frase[a++] = capsula->Cabecalho->nome[j];
                j++;
            }
            frase[a] = '\0';

            no* passagem = capsula;
            erro = ValidaUnicidadeCabecalho(passagem, frase, i, log);

            if(erro)
            {
                if(capsula->Cabecalho->nome[strlen(frase) + 1] == '\0')
                    break;
                if(i+1 != 10)
                    erro=0;
                i++;
            }
            else
            {
                int newc = strlen(frase)-1;
                int troca = 0;

                while(frase[newc] == '_'){
                    newc--;
                    troca++;
                }

                if(troca)
                    frase[newc] = '\0';
                break;
            }
        }

        if(erro)
        {
            if(geralog){
                log = fopen(caminho_arquivo_log, modo_escrita_continua);
                fprintf(log, "Nome de cabecalho repetido: %s\n", capsula->Cabecalho->nome);
                fprintf(log, "Renomeando todos os cabecalhos...\n\n");
                fclose(log);
            }

            Renomear(F);
            break;
        }

        if(i>maior) maior = i;
        capsula = capsula->prox;
    }

    capsula = F->inicio;
    frase[0] = '\0';
    while(capsula != NULL)
    {
        i = 0;
        j = 0;

        while(j<maior && capsula->Cabecalho->nome[i]!='\0')
        {
            if(capsula->Cabecalho->nome[i] == ' ' ||
	       capsula->Cabecalho->nome[i] == '|' ||
	       capsula->Cabecalho->nome[i] == ',' ||
	       capsula->Cabecalho->nome[i] == ';')
            {
                j++;
                if(j<maior)frase[i++] = '_';
            }
            else{	      
	      	frase[i] = capsula->Cabecalho->nome[i];
                i++;
            }
        }
        frase[i] = '\0';

        strcpy(capsula->Cabecalho->nome, frase);

        capsula = capsula->prox;
    }

    return 1;
}

void Renomear(Fila* F)
{
    no* capsula = F->inicio;
    int i=1;
    char temp[100];
    while(capsula != NULL)
    {
        snprintf(temp, 900, ">CDS_%06d", i);
        strcpy(capsula->Cabecalho->nome, temp);
        capsula = capsula->prox;
        i++;
    }
}

int ValidaSequencia(sequencia* seq)
{
    int vetor[26], i;
    memset(vetor, 0, sizeof(int)*26);

    for(i=0; i < seq->tamanho; i++)
    {
        if(seq->sequencia[i]-65 > 25)
            return 0;
        vetor[seq->sequencia[i]-65]++;
    }

    i=26;
    int cont=0;
    while(i--)
    {
        if(vetor[i])
        {
            cont++;
        }
    }
    if(cont == 4 || cont == 5)
    {
        if(vetor[0] && vetor[2] && vetor[6] && (vetor[19] || vetor[20]))
            return 2;
    }

    return 1;
}

int ValidaSequenciaAminoacido(Fila* F, FILE* log){
    if(geralog)
    {
        log = fopen(caminho_arquivo_log, modo_escrita_continua);
        fprintf(log, "Iniciando validação das Sequencias.\n\n");
        fclose(log);
    }

    int resultado=0;
    no* capsula = F->inicio;
    if(capsula == NULL) return 0;
    while(capsula != NULL)
    {
        sequencia* s = capsula->Cabecalho->aminoacidos;

        resultado = ValidaSequencia(s);

        switch(resultado)
        {
            case 2:
                capsula->correto = 0;
                break;
            case 0:
                capsula->correto = 0;
                break;
        }

        capsula = capsula->prox;
    }
    return 1;
}
