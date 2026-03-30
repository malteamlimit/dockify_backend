Use `environment.yml` as the source of truth for runtime dependencies.
The setup is conda-only (no separate `requirements.txt`/pip layer).

Create the environment:

```bash
conda env create -f environment.yml
conda activate dockify
```

Then run api by:
```bash
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```

## Docker (base + app split)

To avoid re-downloading heavy dependencies like `pyrosetta` on every rebuild,
this project uses two Dockerfiles:

- `Dockerfile.base`: heavy and rarely changing dependencies
- `Dockerfile`: lightweight app image with your code

Dependency split:

- `environment.base.yml`: `pyrosetta`, `rdkit`, `openbabel`
- `environment.app.yml`: `fastapi`, `uvicorn`, `sqlmodel`, `python-multipart`

Build base image once (or only when dependency set changes):

```bash
docker build --platform linux/amd64 -f Dockerfile.base -t dockify-base:py310 .
```

Build app image during development:

```bash
docker build --platform linux/amd64 -t dockify-backend:dev .
```

Run app image:

```bash
docker run --platform linux/amd64 --rm -p 8000:8000 dockify-backend:dev
```

Run with live code mount (no rebuild for Python code changes):

```bash
docker run --platform linux/amd64 --rm -it \
  -p 8000:8000 \
  -v "$PWD/app:/app/app" \
  -v "$PWD/input:/app/input:ro" \
  dockify-backend:dev \
  uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```
