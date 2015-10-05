### Reactions file specification

#### Specifica 

Un file *reactions* contiene l'elenco delle reazioni
chimiche associate al sistema.

Il file testuale inizia con due righe che specificano
il numero di reazioni chimiche ed il numero di specie
presenti nel sistema.


    reactions: <numero di reazioni>
    species: <numero di specie>

dopo una riga vuota, seguono `r = <numero di reazioni>`
righe, una per ogni reazione, che indicano quali sono
i reagenti ed i prodotti della reazione.

    <r_1> <r_2> ... <r_s> -> <p_1> <p_2> ... <p_s>

Gli `s = <numero di specie>` reagenti e gli `s = <numero di specie>`
prodotti sono separati da un `->`.

#### Esempio

Un esempio di file *reactions* completo:

    reactions: 2
    species: 3
    
    0 1 1 -> 1 0 0
    1 0 1 -> 1 1 0

#### Limitazioni

- il numero di linee di reazione deve essere uguale ad `r`.
- ogni linea di reazione deve contenere `s` reagenti, un
simbolo `->`, ed `s` prodotti.
- ogni `r_n` può  valere `0` oppure `1`.
- la somma degli `r_n` e la somma dei `p_n` devono valere
`0` oppure `1` oppure `2` (ovvero, in ogni reazione, sono
permessi zero, uno o due reagenti e zero, uno o due prodotti).



