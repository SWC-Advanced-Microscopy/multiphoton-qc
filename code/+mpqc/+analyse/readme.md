# This folder contains various analysis functions

## Converting images to photons
Start with `mpqc.analyse.get_quantalsize_quantalsize_from_file` and `.compute_quantalsize`

```
OUT = mpqc.analyse.get_quantal_size_from_file(fname);
mpqc.analyse.plotPhotonFit(OUT)
```

The above is doing:
```
[im,metadata]=mpqc.tools.scanImage_stackLoad(fname,false); %Do not subtract the offset
```

now fits the model:
```
[result,dataForFit] = mpqc.analyse.compute_quantalsize(im);
```

and check the fit:
```
mpqc.analyse.plotPhotonFit(result)
```
