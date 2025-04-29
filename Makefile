ENGINE = podman
TARGET = ngm

.PHONY: help local deploy clean build_container run_container

help: # show help for each of the Makefile recipes
	@grep -E '^[a-zA-Z0-9 _-]+:.*#'  Makefile | while read -r l; do printf "\033[1;32m$$(echo $$l | cut -f 1 -d':')\033[00m:$$(echo $$l | cut -f 2- -d'#')\n"; done

local: # run app in local environment
	streamlit run app.py

build_container: # build container locally
	$(ENGINE) build -t $(TARGET) -f Dockerfile

run_container: # run container locally
	$(ENGINE) run -p 8501:8501 --rm $(TARGET)

deploy: manifest.json requirements.txt # deploy via Posit Connect
	rsconnect deploy \
		manifest manifest.json \
		--title $(TARGET)

manifest.json requirements.txt: app.py ngm/app.py ngm/linalg.py ngm/__init__.py pyproject.toml poetry.lock
	rm -f requirements.txt
	rsconnect write-manifest streamlit . \
		--overwrite \
		-x ".github/**" \
		-x .gitignore \
		-x .pre-commit-config.yaml \
		-x ".pytest_cache/**" \
		-x ".ruff_cache/**" \
		-x .secrets.baseline \
		-x ".vscode/**" \
		-x Dockerfile \
		-x LICENSE \
		-x Makefile \
		-x README.md \
		-x "docs/**" \
		-x mkdocs.yml \
		-x "ngm/__pycache__/**" \
		-x poetry.lock \
		-x pyproject.toml \
		-x "tests/**"

clean: # remove Posit Connect intermediates
	rm -f manifest.json requirements.txt
