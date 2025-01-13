ENGINE = podman
TARGET = ngm

.PHONY: run build_container run_container

run:
	streamlit run ngm/app.py

build_container:
	$(ENGINE) build -t $(TARGET) -f Dockerfile

run_container:
	$(ENGINE) run -p 8501:8501 --rm $(TARGET)
