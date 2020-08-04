typedef struct fila Fila;

typedef struct sequencia {char *sequencia;int tamanho;} sequencia;
typedef struct cabecalho {char *nome;struct sequencia *aminoacidos;} cabecalho;

typedef struct bloco {cabecalho *Cabecalho;struct bloco *prox; int correto;} no;
int geralog;
int geraArqSaida;
int gerarNomeLog;
int gerarNomeSaida;

void SetaValores(char* c, char* a, char* b);

void Esvazia(Fila *F);

void Destroi(Fila *F);

Fila* Cria();

void Entra(Fila *F, cabecalho *Cab, sequencia *Seq);

void Sai(Fila *F, cabecalho *Cab, sequencia *Seq);

void Inicio(Fila* F, cabecalho *Cab);

int Vazia(Fila *F);

int Tamanho(Fila *F);

int DescarregaFilaEmArquivo(Fila* F, FILE* arquivo, FILE* log);

int DescarregaFilaEmTela(Fila* F, FILE* log);

int ValidaCabecalhoAminoacido(Fila* F, FILE* log);

int ValidaSequenciaAminoacido(Fila* F, FILE* log);
