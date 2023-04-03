FROM python:3.9.13
COPY requirements.txt requirements.txt
RUN pip install -r requirements.txt
#COPY data/ /data
COPY van.py /tmp/
CMD ["python", "/tmp/van.py"]
