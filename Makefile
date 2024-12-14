ifndef ENGINE
ENGINE = podman
endif

ifndef TARGET
TARGET = ngm
endif

.PHONY: streamlit build_container

streamlit:
	streamlit run scripts/widget.py

build_container:
	$(ENGINE) build -t $(TARGET) -f Dockerfile

run_container:
	$(ENGINE) run -p 8501:8501 --rm $(TARGET)
