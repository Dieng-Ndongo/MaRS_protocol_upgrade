import os
import sys
import subprocess
from pathlib import Path


ENV_NAME = "pipeline_env"
ENV_FILE = "environment.yml"
MINICONDA_DIR = os.path.expanduser("~/pipeline/miniconda3")

# Contenu par défaut de l'environnement
DEFAULT_ENV = f"""name: {ENV_NAME}
channels:
  - bioconda
  - conda-forge
  - anaconda
  - conda
  - gstr
  - defaults
dependencies:
# Python
  - python=3.10
  - pip
  - pip:
    - streamlit
# Bioinformatique / outils NGS
  - fastqc
  - multiqc
  - bwa
  - samtools
  - bcftools
  - bedtools
  - gatk4
  - trimmomatic
  - cutadapt
  - sickle
  - seqtk
  - pysam
  - bbmap # pour bbduck
  - freebayes
  - vardict
  - vardict-java
  - perl
  - snpeff 
  - snpsift
  - picard
# Java pour snpEff
  - openjdk=21
# Librairies Python utiles
  - pandas
  - numpy 
  - pyfaidx
  - tqdm
  - plotly
  - matplotlib
  - seaborn
  - pysam
  - biopython
  - reportlab
# Variables d'environnement pour Java
variables:
   JAVA_HOME: $CONDA_PREFIX
"""

def run_cmd(cmd, capture_output=False, **kwargs):
    """Exécute une commande shell et affiche la commande exécutée"""
    print(f"[CMD] {' '.join(cmd)}")

    # Si tu veux capturer la sortie (logs internes)
    if capture_output:
        return subprocess.run(cmd, check=True, text=True, capture_output=True, **kwargs)
    else:
        # Affiche en direct les sorties du processus
        return subprocess.run(cmd, check=True, text=True, stdout=sys.stdout, stderr=sys.stderr, **kwargs)

def check_conda():
    """Vérifie si conda est installé"""
    try:
        run_cmd(["conda", "--version"])
        return True
    except Exception:
        return False
    
def install_miniconda():
    """Télécharge et installe Miniconda localement"""
    if not os.path.exists(MINICONDA_DIR):
        print("[INFO] Téléchargement et installation de Miniconda...")
        run_cmd([
            "wget", "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh",
            "-O", "miniconda.sh"
        ])
        run_cmd(["bash", "miniconda.sh", "-b", "-p", MINICONDA_DIR])
    else:
        print("[INFO] Miniconda déjà installé.")

    os.environ["PATH"] = f"{MINICONDA_DIR}/bin:" + os.environ["PATH"]

run_cmd(["conda", "config", "--add", "channels", "conda-forge"])
run_cmd(["conda", "config", "--add", "channels", "bioconda"])
run_cmd(["conda", "config", "--set", "channel_priority", "strict"])


def create_or_update_env():
    """Crée ou met à jour l'environnement conda depuis environment.yml"""
    if not Path(ENV_FILE).exists():
        print(f"[WARN] Fichier {ENV_FILE} introuvable, création d’un fichier par défaut...")
        with open(ENV_FILE, "w", encoding="utf-8") as f:
            f.write(DEFAULT_ENV.strip() + "\n")

    # Vérifie si l'environnement existe déjà
    result = subprocess.run(["conda", "env", "list"], capture_output=True, text=True)
    env_exists = any(ENV_NAME in line for line in result.stdout.splitlines())

    if env_exists:
        print(f"[INFO] Mise à jour de l'environnement {ENV_NAME}...")
        run_cmd(["conda", "env", "update", "-f", ENV_FILE, "--prune"])
    else:
        print(f"[INFO] Création de l'environnement {ENV_NAME}...")
        run_cmd(["conda", "env", "create", "-f", ENV_FILE])


# def run_pipeline():
 #   """Exécute le pipeline dans l'environnement conda"""
 #   print("[INFO] Lancement du pipeline dans l'environnement conda...")
 #   run_cmd(["conda", "run", "-n", ENV_NAME, "python", "pipeline_python.py"],capture_output=False )

def run_cmd_live(cmd, env=None):
    """Exécute une commande shell et affiche la sortie en direct"""
    print(f"[CMD] {' '.join(cmd)}")
    process = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env  # permet de passer un environnement modifié
    )
    for line in process.stdout:
        print(line, end="")  # affiche chaque ligne en direct
    process.wait()
    if process.returncode != 0:
        raise subprocess.CalledProcessError(process.returncode, cmd)


def main():
    if not check_conda():
        print("[INFO] Conda non trouvé. Installation de Miniconda...")
        install_miniconda()  
    create_or_update_env()

    # Trouver le chemin de base de conda
    conda_prefix = subprocess.run(
        ["conda", "info", "--base"], capture_output=True, text=True
    ).stdout.strip()

    if not conda_prefix:
        print("[ERROR] Conda n'est pas installé ou non trouvé.")
        sys.exit(1)

    # Chemin vers l'environnement conda cible
    env_path = os.path.join(conda_prefix, "envs", ENV_NAME)
    env = os.environ.copy()
    env["PATH"] = os.path.join(env_path, "bin") + os.pathsep + env["PATH"]

    # Lancer le pipeline dans l'environnement conda
    cmd = ["python", "pipeline_python.py"]
    run_cmd_live(cmd, env=env)

if __name__ == "__main__":
    main()
