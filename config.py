from pathlib import Path
from sqlalchemy import create_engine

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"

MOLECULES_FILE = DATA_DIR / "molecules.tsv"
ADVERSE_EFFECTS_FILE = DATA_DIR / "adverseEffects.tsv"

DATABASE_NAME = "drug_database.db"
DATABASE_PATH = BASE_DIR / DATABASE_NAME
DATABASE_URL = f"sqlite:///{DATABASE_PATH}"

ENGINE = create_engine(DATABASE_URL)
API_URL = "http://127.0.0.1:8000"