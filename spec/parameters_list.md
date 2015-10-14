### CUDA functions parameter list

#### Generale

```
function(
    int * state,              // state array (len = sbc * spc)
    int * reactants,          // reactants array (len = rc * spc)
    int * products,           // products array  (len = rc * spc)
    int * topology,           // topology array  (len = 6 * sbc)
    int sbc,                  // subvolumes count
    int spc,                  // species count
    int rc,                   // reactions count
    int sbi,                  // subvolume index
    int spi,                  // specie index
    int ri,                   // reaction index
    float * rate_matrix,      // rate matrix     (len = 3 * sbc)
    float * rrc,              // reactions rates constants  (len = rc)
    float * drc,              // diffusion rates constants  (len = spc)
)
```

