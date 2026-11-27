"""Совместимое имя проверки; актуальный протокол первого этапа не перезаписывается."""
from pathlib import Path
from runpy import run_path
run_path(str(Path(__file__).with_name('verify_stage1_audit.py')), run_name='__main__')
