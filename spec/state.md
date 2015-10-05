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
- il file deve concludersi con al piu' 1 *whitespace*, che
segue le linee di reazione.



