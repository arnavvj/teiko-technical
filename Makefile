VENV   := .venv
PYTHON := $(VENV)/bin/python

.PHONY: setup pipeline dashboard clean

setup:
	python3 -m venv $(VENV)
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) analysis.py

dashboard:
	$(PYTHON) -m streamlit run dashboard.py --server.headless true --server.port 8501

clean:
	rm -rf cell_counts.db outputs
