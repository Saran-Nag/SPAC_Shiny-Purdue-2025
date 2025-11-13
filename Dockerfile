FROM spac-shiny-base:latest

WORKDIR /app

COPY . .

EXPOSE 8000

CMD ["sh", "-c", "python -u -m shiny run app.py --host 0.0.0.0 --port 8000 2>&1 | sed -u 's|0\\.0\\.0\\.0|localhost|g'"]