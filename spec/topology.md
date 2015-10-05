### State file specification

#### Specifica 

Un file *topology* contiene una descrizione della topologia
del sistema.

Il file testuale inizia con una riga che specifica
il numero di sottovolumi di cui la topologia del sistema
è composta.

    subvolumes: <numero di sottovolumi>

dopo una riga vuota, seguono `N = <numero di sottovolumi>`
righe, una per ogni sottovolume, che indicano quali 
sottovolumi sono adiacenti al sottovolume corrente.

    <n>: <s_1> <s_2> ... <s_x>

`n` è il numero del sottovolume, e deve essere compreso tra
`0` e `N-1`. `s_x` indica che il sottovolume numero `x` è
adiacente al sottovolume `n`. 

#### Esempio

Un esempio di file *topology* completo:

    subvolumes: 4
    
	0: 1
	1: 0 2
	2: 1 3
	3: 2

Che corrisponde a

      +-----+-----+-----+-----+
     /     /     /     /     /|
    +-----+-----+-----+-----+ |
    |  0  |  1  |  2  |  3  | / 
    +-----+-----+-----+-----+-

#### Limitazioni

- il numero di linee di sottovolume deve essere uguale ad `N`.
- ogni linea di sottovolume deve contenere al più `6` valori
compresi tra `0` ed `N-1`.




