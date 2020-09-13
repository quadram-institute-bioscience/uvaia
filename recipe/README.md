## Integration recipes  

The instructions below are not useful for the final user. 
bioconda needs `meta.yaml` and `build.sh` on a separate folder.
```bash
conda-build super_distance
# after running bootstrap on bioconda-recipes
bioconda-utils lint recipes super_distance/meta.yaml --git-range master
bioconda-utils build --packages super_distance --testonly recipes/ config.yml
bioconda-utils build recipes config.yml --docker --mulled-test --git-range master
```
Docker is set up remotely, with `Dockerfile` file on root of repo (and uses my personal fork github.com/leomrtns)
```bash
docker pull leomrtns/super_distance
```
travis.ci depends on .travis.yml at repo root. 
