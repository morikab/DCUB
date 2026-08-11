import sys
from pathlib import Path

# `modules` (app/modules/) is imported everywhere in this codebase as a
# top-level package — `from modules import models`, not `from app.modules
# import models` — matching how api_server.py and the analysis notebooks
# already run with app/ on sys.path (see this repo's CLAUDE.md: "Notebooks
# add the project root to sys.path"). Do the same here so tests can import
# the real backend code without duplicating it.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
