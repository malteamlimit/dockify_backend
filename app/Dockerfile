FROM --platform=linux/amd64 python:3.9

WORKDIR /code

COPY ./requirements.txt /code/requirements.txt

# Install packages from PyPI directly
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r /code/requirements.txt

# Copy application code
COPY . /code

# Command to run the application
CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "80"]