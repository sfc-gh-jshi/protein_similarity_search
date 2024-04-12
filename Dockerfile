FROM python:3.10
WORKDIR /app
ADD ./requirements.txt /app/
RUN pip install --no-cache-dir -r requirements.txt
ADD ./ /app
EXPOSE 8080
CMD ["streamlit", "run", "--logger.level=debug",  "--server.port=8080",  "--server.address=0.0.0.0",  "--server.runOnSave=true",  "--server.fileWatcherType=poll", "protein_streamlit.py"]