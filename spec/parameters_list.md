### CUDA functions parameter list

#### Generale

```
function(
    int * state,              // state array (len = SBC * SPC)
    int * reactants,          // reactants array (len = RC * SPC)
    int * products,           // products array  (len = RC * SPC)
    int * topology,           // topology array  (len = 6 * SBC)
    int sbi,                  // subvolume index
    int spi,                  // specie index
    int ri,                   // reaction index
    float * rate_matrix,      // rate matrix     (len = 3 * SBC)
    float * rrc,              // reactions rates constants  (len = RC)
    float * drc,              // diffusion rates constants  (len = SPC)
)
```

