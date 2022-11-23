FROM python

COPY ./ /app
WORKDIR /app
RUN ls

RUN pip install --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt

CMD ["python", "main.py"]
