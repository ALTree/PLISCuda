### Rate matrix representation 

La parte del sistema concernente tutto ciò che riguarda
la gestione dei *reaction rates* è composta da tre
componenti: una *rate matrix*, che contiene i tassi di
reazione e diffusione di ogni sottovolume, (e la loro
somma) e due array *reaction rates constants* e
*diffusion rates constants* che contengono le costanti
di reazione e diffusione specifiche di ogni specie
contenuta nel sistema.

#### Specifica 

`{reaction, diffusion}_rates_costants` sono, per il momento,
*hard-coded* nel sorgente, all'interno di due array

    []float = {r_0, r_1, ... r_x}
    []float = {s_0, s_1, ... s_x}

dove `x = <numero di sottovolumi> - 1`

**TODO:** lettura di `{reaction, diffusion}_rates_constants` da file(?)

**TODO:** *rates constants* differenti per ogni sottovolume(?)

La matrice dei tassi di reazione contiene la somma dei
*reaction rates* (`R`) di tutte le reazioni, la somma di tutti
i *diffusion rates* (`S`), e la somma delle due (`R + S`), e
viene calcolata a `run-time` a partire dalle *reaction* e *diffusion*
*rates constants* e dal numero di molecole presenti nel sottovolume.

Il *reaction rate* di una reazione `r` in un sottovolume
`s` dipende dal tipo di reazione e dalla quantità e dal tipo di molecole
presenti nel sottovolume. In particolare:

1. Per le reazioni di tipo `uni`, il *reaction rate* è uguale al numero
di molecole della specie coinvolta moltiplicato per la *reaction rate constant*
associata alla reazione
2. Per le reazioni di tipo `bi_same`, posto `N` il numero di molecole della
specie coinvolta presenti nel sottovolume, il *reaction rate* è uguale a
`0.5 * (N) * (N-1) * reaction_rate_constant`
3. Per le reazioni di tipo `bi_diff`, posti `N1` ed `N2` i numeri di molecole
della prima e della seconda specie coinvolte, il *reaction rate* èè uguale a
`N1 * N2 * reaction_rate_constant`


#### Esempio

    float r_rates_constants[] = {0.1, 0.1, 0.2}
    float d_rates_constants[] = {0.5, 0.5, 0.6}

    | R | S | R+S |
    |   |   |     |
    |   |   |     |
    |...|...| ... |

#### Limitazioni

- *reaction rates constants* e *diffusion rates constants*
devono essere nulli o positivi.

**TODO:** strettamente positivi(?)

#### Rappresentazione in memoria

*reaction rates constants* e *diffusion rates constants* sono memorizzati in due
zone di memoria contigue (due array mono-dimensionali).

La *rate matrix* è memorizzata in una zona di memoria contigua
(i.e. un array mono-dimensionale) di `3n` (con `n = <numero di sottovolumi>`)
contenente nelle celle dalla `0` alla `n-1` tutti le *reaction rates constants*,
nelle celle dalla `n` alla `2n-1` le *diffusion rates costants*, e nelle
celle dalla `2n` alla `3n-1` le somme delle due.

