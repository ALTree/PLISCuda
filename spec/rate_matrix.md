### Rate matrix representation 

La parte del sistema concernente tutto ciò che riguarda
la gestione dei *reaction rates* è composta da tre
componenti: una *rate matrix*, che contiene i tassi di
reazione e diffusione di ogni sottovolume, (e la loro
somma) e due array *reaction rates* e *diffusion rates*
che contengono le costanti di reazione e diffusione
specifiche di ogni specie contenuta nel sistema.

#### Specifica 

`{reaction, diffusion}_rates` sono, per il momento,
*hard-coded* nel sorgente, all'interno di due array

    []float = {r_0, r_1, ... r_x}
    []float = {s_0, s_1, ... s_x}

dove `x = <numero di sottovolumi> - 1`

**TODO:** lettura di `{reaction, diffusion}_rates` da file(?)
**TODO:** *rates* differenti per ogni sottovolume(?)

La matrice dei tassi di reazione contiene la somma dei
*reaction rates* (`R`) di tutte le reazioni, la somma di tutti
i *diffusion rates* (`S`), e la somma delle due (`R + S`), e
viene calcolata a `run-time` a partire dai *reaction* e *diffusion*
*rates* e dal numero di molecole presenti nel sottovolume.

#### Esempio

    float r_rates[] = {0.1, 0.1, 0.2}
    float d_rates[] = {0.5, 0.5, 0.6}

    | R | S | R+S |
    |   |   |     |
    |   |   |     |
    |...|...| ... |

#### Limitazioni

*reaction rates* e *diffusion rates* devono essere nulli
o positivi.

*TODO:* strettamente positivi(?)

#### Rappresentazione in memoria

*reaction rates* e *diffusion rates* sono memorizzati in due
zone di memoria contigue (due array mono-dimensionali).

La *rate matrix* è memorizzata in una zona di memoria contigua
(i.e. un array mono-dimensionale) di `3n` (con `n = <numero di sottovolumi>`)
contenente nelle celle dalla `0` alla `n-1` tutti i *reaction rates*,
nelle celle dalla `n` alla `2n-1` i *diffusion rates*, e nelle
celle dalla `2n` alla `3n-1` le somme dei due.

