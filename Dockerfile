FROM spac-shiny-base:latest

WORKDIR /app

COPY . .

EXPOSE 8000

CMD ["python", "-m", "shiny", "run", "app.py", "--host", "0.0.0.0", "--port", "8000"]