### State file specification

#### Specifica 

Un file *state* contiene lo stato iniziale del sistema.

Il file testuale inizia con due righe che specificano
il numero di sottovolumi ed il numero di specie presenti
nel sistema.


    subvolumes: <numero di sottovolumi>
    species: <numero di specie>

dopo una riga vuota, seguono `N = <numero di sottovolumi>`
righe, una per ogni sottovolume, che indicano quali specie
chimiche sono presenti in quel sottovolume.

    <n>: <s_1> <s_2> ... <s_s>

`n` è il numero del sottovolume, e deve essere compreso tra
`0` e `N-1`. `s_x` indica la quantità di molecole della specie
`x` presenti nel sottovolume numero `n`. Gli `s_x` devono essere
in numero pari ad `s = <numero di specie>` indicato nell'header.

#### Esempio

Un esempio di file *state* completo:

    subvolumes: 3
    species: 4
    
    0: 0 0 100
    1: 0 100 0
	2: 100 0 0

#### Limitazioni

- il numero di linee di sottovolume deve essere uguale ad `N`.
- ogni linea di sottovolume deve contenere `s` valori interi.
- ogni `s_n` deve essere un numero intero positivo.


#### Rappresentazione in memoria

Lo stato del sistema è memorizzato in una zona di memoria
contigua (i.e. un array mono-dimensionale).

Dato un sistema avente `n` sottovolumi ed `s` specie, l'array dello
stato contiene il numero di molecole presenti in ogni sottovolume,
ordinato per *specie*. Le celle di memoria dalla `0` alla `n-1` contengono
le `n` quantità presenti della specie `0`; le celle di memoria
dalla `n` alla `2n-1` contengono le `n` quantità presenti della
specie `1`, e così via.

Ad esempio, il sistema a `3` sottovolumi avente stato

    0: 0 0 5
    1: 0 5 0
	2: 5 0 4

sarà rappresentato in memoria come

    []state = {0 0 5 0 5 0 5 0 4}


