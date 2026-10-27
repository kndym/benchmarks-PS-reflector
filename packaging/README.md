# Packaging and build environments

[docker/](docker/) contains the repository's Docker build environment.
Run builds from the repository root so the Docker context includes the source:

```bash
docker build -f packaging/docker/Dockerfile .
```

Packaging files belong here; computed arrays, plots, and logs belong in the
appropriate [results category](../results/). Older write-up ZIP exports are
preserved under [archive](../archive/writeup_exports_oct_27_26/).

<!-- directory-index:start -->

## Directory guide

[Parent directory](../README.md)

| Folder | Contents |
| --- | --- |
| [docker/](docker/README.md) | Supporting files grouped by purpose. Follow the directory links below and the parent README for context. |

<!-- directory-index:end -->
