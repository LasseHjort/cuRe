### cuRe: An R package for estimating cure models and loss of lifetime ###


Installation
-----------

The following packages needs to be installed

```
pkgs <- c("numDeriv", "relsurv")
install.packages(pkgs)
```

The development branch of the rstpm2 package has to be installed in order to have proper generation of initial values.
```
library(githubinstall)
gh_install_packages("rstpm2", ref = "develop")
```

The cuRe package is then installed by
```
gh_install_packages("cuRe")
```


