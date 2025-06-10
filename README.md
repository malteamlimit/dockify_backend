pyrosetta must be installed in the same environment as this project.

Start by installing the requirements:

```bash
pip install -r requirements.txt
```

Then run api by:
```bash
uvicorn app.main:app --reload --host 0.0.0.0 --port 8000
```