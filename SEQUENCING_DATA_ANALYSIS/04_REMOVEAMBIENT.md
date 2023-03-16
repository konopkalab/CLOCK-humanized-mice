# REMOVE AMBIENT RNA USING CELLBENDER
- python/3.7.x-anaconda
- cellbender

## P07 | WT
```{sh}
ls | grep -P "^YL_P07_WT" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000'
```

## P07 | KO
```{sh}
ls | grep -P "^YL_P07_KO" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000'
```

## P07 | HU
```{sh}
ls | grep -P "^YL_P07_HU" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000'
```
## P56 | WT
```{sh}
ls | grep -P "^YL_P56_WT" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000'
```

## P56 | KO
```{sh}
ls | grep -P "^YL_P56_KO" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000'
```

## P56 | HU
```{sh}
ls | grep -P "^YL_P56_HU" | xargs -I % -n 1 -P 48 sh -c 'echo %; cellbender remove-background --input % --output %/CellBender_Out.h5 --expected-cells 10000'
```

-----
