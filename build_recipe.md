# To properly build the Docker images, run these in sequence:

```bash
$ docker build -t ninai/pipline:base -f Dockerfile.base .
$ docker build -t ninai/pipeline .
$ docker build -t ninai/pipeline:tf -f Dockerfile.tf-gpu .
```

