
# Makefile for installing the CRISPR spacer identification tool

.PHONY: install setup_scripts

# Path to the BLAST binaries
BLAST_PATH=/usr/bin

# Directory to place the scripts
SCRIPT_DIR=/usr/local/bin

install: setup_scripts
	@echo "Installing Python dependencies..."
	pip install -r requirements.txt
	@echo "Checking if BLAST is installed..."
	@if [ -x "$(BLAST_PATH)/blastn" ]; then 		echo "BLAST is already installed."; 	else 		echo "BLAST is not installed. Please install BLAST from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/"; 	fi
	@echo "Installation completed."

setup_scripts:
	@echo "Setting up the scripts..."
	# Copy and set up the Python scripts
	cp createBlastDB.py $(SCRIPT_DIR)
	cp spacePhage.py $(SCRIPT_DIR)
	chmod +x $(SCRIPT_DIR)/createBlastDB.py
	chmod +x $(SCRIPT_DIR)/spacePhage.py
	# Copy and set up the Bash script
	cp spacePhage.sh $(SCRIPT_DIR)/spacePhage
	chmod +x $(SCRIPT_DIR)/spacePhage
	@echo "Scripts have been set up."
