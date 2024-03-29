### Rate matrix representation 

La parte del sistema concernente tutto ciò che riguarda
la gestione dei *reaction rates* è composta da cinque
componenti: una *rate matrix*, che contiene i tassi di
reazione e diffusione di ogni sottovolume, (e la loro
somma), due array *reaction rates constants* e
*diffusion rates constants* che contengono le costanti
di reazione e diffusione specifiche di ogni reazione e
di ogni specie considerate nel sistema, e due array
*react_rates_array* e *diff_rates_array* che contengono
i tassi di reazione e diffusione (rispettivamente per ogni
reazione ed ogni specie) in ogni sottovolume del sistema.

#### Specifica 

`{reaction, diffusion}_rates_costants` sono letti da file e
memorizzati all'interno di due array

    []float = {r_0, r_1, ... r_x}
    []float = {s_0, s_1, ... s_x}

dove `x = <numero di sottovolumi> - 1`

**TODO:** *rates constants* differenti per ogni sottovolume(?)

`{react, diff}_rates_array` sono calcolati a run time e memorizzati
all'interno di due array. `react_rates_array` ha dimensioni pari
al numero di sottovolumi moltiplicato per il numero di reazioni,
mentre `diff_rates_array` ha dimensioni pari al numero di sottovolumi
moltiplicato per il numero di specie.

La *rate matrix* è composta da 3 colonne e un numero di righe pari
al numero di sottovolumi. Per ogni riga (e quindi per ogni sottovolume),
sono memorizzate:

1. La somma dei *reaction rates* di tutte le reazioni (`R`)
2. La somma dei *diffusion rates* di tutte le specie (`S`)
3. `R + S`

La tabella viene inizialmente computata e successivamente aggiornata
a `run-time`.

##### Calcolo reaction rates

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

La somma su tutte le reazioni dei *reaction rates* (che chiamiamo `R`), viene
memorizzata nella prima colonna della *rate matrix*.

##### Calcolo diffusion rates

Il *diffusion rate* di una specie fornisce una stima numerica della propensità
che tale molecola ha di saltare in uno dei sottovolumi adiacenti a quello in
cui si trova (ovvero, di muoversi).

Il *diffusion rate* di una singola specie è dato dal prodotto della
*diffusion rate constant* per quella specie moltiplicato per il numero
di molecole di quella specie presenti nel sottovolume.

La somma su tutte le specie dei *diffusion rates* moltiplicata per il numero
di sottovolumi adiacenti a quello considerato, (che chiamiamo `S`), viene
memorizzata nella seconda colonna della *rate matrix*.


#### Esempio

Dato un sistema con 3 reazioni e 4 specie e 3 sottovolumi:

    float r_rates_constants[] = {0.1, 0.1, 0.2}
    float d_rates_constants[] = {0.5, 0.5, 0.5, 0.6}

       | R | S | R+S |
    s0 | . | . |  .  |
    s1 | . | . |  .  |
    s2 | . | . |  .  |

I rate di reazioni sono in numero pari al numero di reazioni, i
rate di diffusione sono in numero pari a quello delle specie,
la *rate matrix* contiene una riga per ogni sottovolume di cui il
sistema è composto.

#### Limitazioni

- *reaction rates constants* e *diffusion rates constants*
devono essere nulli o positivi.

**TODO:** strettamente positivi(?)

#### Rappresentazione in memoria

*reaction rates constants* e *diffusion rates constants* sono memorizzati in due
zone di memoria contigue (due array mono-dimensionali).

`react_rates_array` e `diff_rates_array` sono organizzati rispettivamente per
reazione e per specie.

I *react rates* della reazione `0` occupano le prime `sbc`
posizioni nell'array (in ordine di sottovolume), i *rates* della reazione `1`
occupano le posizioni dalla `sbc` alla `2*sbc-1`, e così via.

I *diff rates* della specie `0` occupano le prime `sbc`
posizioni nell'array (in ordine di sottovolume), i *rates* della specie `1`
occupano le posizioni dalla `sbc` alla `2*sbc-1`, e così via.

La *rate matrix* è memorizzata in una zona di memoria contigua
(i.e. un array mono-dimensionale) di `3n` (con `n = <numero di sottovolumi>`)
contenente nelle celle dalla `0` alla `n-1` tutti le *reaction rates constants*,
nelle celle dalla `n` alla `2n-1` le *diffusion rates costants*, e nelle
celle dalla `2n` alla `3n-1` le somme delle due.

