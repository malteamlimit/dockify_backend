FROM dockify-base:py310

WORKDIR /app

COPY environment.app.yml ./environment.app.yml
RUN conda env update -n dockify -f environment.app.yml && conda clean -afy

COPY app ./app
COPY input ./input

RUN mkdir -p data app/static/poses app/static/previews

EXPOSE 8000

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
