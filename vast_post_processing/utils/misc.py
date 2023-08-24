"""

 Miscellaneous Utilities for VAST Post Processing

"""
from pathlib import Path
from typing import Tuple

def get_epoch_directory(image_path: Path) -> str:
    return image_path.parent.name

def get_field_and_sbid(image_path: Path) -> Tuple[str,int]:
    _, _, field, sbid_str, *_ = image_path.name.split(".")
    sbid = int(sbid_str[2:])
    return field, sbid