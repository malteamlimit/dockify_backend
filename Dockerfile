FROM continuumio/miniconda3

RUN conda create -n dockify python=3.12 -y

WORKDIR /app

COPY environment.yml .

RUN conda env update -n dockify -f environment.yml

COPY requirements.txt .
COPY app .

RUN conda run -n dockify pip install --no-cache-dir -r requirements.txt


CMD ["conda", "run", "--no-capture-output", "-n", "dockify", "uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
